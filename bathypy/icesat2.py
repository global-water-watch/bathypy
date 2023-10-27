from glob import glob
import os
from pathlib import Path
import re
from logging import getLogger, Logger
from datetime import datetime, date

import pandas as pd
import geopandas as gpd
import numpy as np


from gwwstorage.open_altimetry import (
    get_track_ids,
    get_overpasses,
    get_icesat2_data_from_OA_api,
)
from gwwstorage.utils import (
    bounding_box_tiles,
    points_in_polygon_parallel,
    buffer_geometry,
)

# functions related to elevation detection and time series

DATE_REGEX: str = r"\d{7}_\d+_(\d{4}-\d{2}-\d{2}).parquet"
FID_REGEX: str = r"(\d{7})_\d+_\d{4}-\d{2}-\d{2}.parquet"

logger: Logger = getLogger(__name__)


def parse_atl03(photon_data: dict, track_id: str, track_date: str) -> list:
    """Parses a dictionary containing data on the atl03 icesat2 product. The dictionary must
    be from the response of the open altimetry API request for requesting icesat2 data.

    :param photon_data: Dictionary containing icesat2 atl03 data
    :type photon_data: dict
    :param track_id: the id of the icesat2 track
    :type track_id: str
    :param track_date: the date of the track
    :type track_date: str
    :return: a list of dict with the icesat2 data
    :rtype: list of dicts
    """
    rows = []

    for beam in photon_data:
        beam_name = beam["beam_name"]
        for s in beam["series"]:
            # every series has name (confidence)
            series_name = s["name"]
            for o in s["data"]:
                # add rows
                row = {
                    "track_id": track_id,
                    "date": track_date,
                    "beam": beam_name,
                    "series": series_name,
                    "lon": round(o[1], 6),
                    "lat": round(o[0], 6),
                    "h": o[2],
                }

                rows.append(row)
    return rows


def parse_atl08(photon_data: dict, track_id: str, track_date: str) -> list:
    """Parses a dictionary containing data on the atl08 icesat2 product. The dictionary must be
    from the response of the open altimetry API request for requesting icesat2 data.

        :param photon_data: Dictionary containing icesat2 atl08 data
        :type photon_data: dict
        :param track_id: the id of the icesat2 track
        :type track_id: str
        :param track_date: the date of the track
        :type track_date: str
        :return: a list of dict with the icesat2 data
        :rtype: list of dicts
    """

    rows = []
    for beam in photon_data["series"]:
        beam_name = beam["beam"]
        isStrongBeam = beam["isStrongBeam"]
        for lle in beam["lat_lon_elev_canopy"]:
            row = {
                "track_id": track_id,
                "date": track_date,
                "beam": beam_name,
                "isStrongBeam": isStrongBeam,
                "lon": round(lle[1], 6),
                "lat": round(lle[0], 6),
                "h": lle[2],
            }
            rows.append(row)
    return rows


def get_bounding_box(aoi: gpd.GeoDataFrame, icesat2_product: str = None) -> list:
    """Creates a bounding box for a given geodataframe. If the icesat product is
    atl03 the size of the bounding box is checked. If the size is bigger than 1
    degree on any axis, the bounding box is cut up in smaller tiles.

    :param aoi: a geodataframe containing a reservoir geometry.
    :type aoi: gpd.GeoDataFrame
    :return: list of bounding boxes, if just one bounding box, the function still
        returns a list of size one.
    """
    minx, miny, maxx, maxy = aoi.loc[0, "geometry"].bounds
    bounding_box = (minx, miny, maxx, maxy)
    if (maxx - minx + maxy - miny) > 1 and icesat2_product == "atl03":
        bounding_box = bounding_box_tiles(bounding_box, step_size=0.5)
    else:
        bounding_box = [bounding_box]
    return bounding_box


def get_clip_feature(gdf: gpd.GeoDataFrame) -> np.ndarray:
    """Creates a slightly larger (200m buffer) simplified geometry from the
    input geometry in the form of a numpy array with coordinates.

    :param gdf: a geodataframe containing a reservoir geometry.
    :type gdf: gpd.GeoDataFrame
    :return: a polygon as an array of points
    :rtype: np.ndarray
    """

    geoseries = gdf.loc[:, "geometry"]
    proj_crs = geoseries.estimate_utm_crs()
    gdf = gdf.to_crs(proj_crs)
    gdf.loc[0, "geometry"] = gdf.loc[0, "geometry"].buffer(200).simplify(50)
    gdf = gdf.to_crs(4326)
    geom_array = np.array(gdf.loc[0, "geometry"].exterior.coords)
    return geom_array


def filter_icesat2_data(
    df: pd.DataFrame, clip_geom: np.array, lon_col: str = "lon", lat_col: str = "lat"
) -> pd.DataFrame:
    """Filters a dataframe containing icesat2 data by calling
    points_in_polygon_parallel function to check if the icesat2 points are within
    a given clip geometry

    :param df: dataframe containing icesat2 points
    :type df: pd.DataFrame
    :param clip_geom: an array of points representing a polygon
    :type clip_geom: np.array
    :param lon_col: name of column with longitude value, defaults to 'lon'
    :type lon_col: str
    :param lat_col: name of column with latitude value, defaults to 'lat'
    :type lat_col: str
    :return: a filtered dataframe with points that are within in the clip geometry
    :rtype: pd.DataFrame
    """
    coords = np.array(list(zip(df[lon_col], df[lat_col]))).astype(np.float32)
    bool_array = points_in_polygon_parallel(coords, clip_geom)
    return df[bool_array]


def process_icesat2_for_gee(df: pd.DataFrame, icesat2_product: str) -> pd.DataFrame:
    """Processes a dataframe containing icesat2 points so that the resulting dataframe
    is ready for uploading to Google Earth Engine

    :param df: dataframe containing icesat2 points
    :type df: pd.DataFrame
    :return: A dataframe with renamed columns and encoded series data
    :rtype: pd.DataFrame
    """
    if icesat2_product == "atl03":
        df.series = df.series.map(
            {"High": 0, "Medium": 1, "Low": 2, "Buffer": 3, "Noise": 4}
        )
    df.date = pd.to_datetime(df.date).values.astype(np.int64) // 10**6

    df = df.rename(
        columns={
            "lon": "longitude",
            "lat": "latitude",
            "date": "system:time_start",
        }
    )
    return df


def drop_redundant_overpasses(overpasses: list, latest_date: str) -> list:
    """Drops the overpasses (track date, track id) from the given list of overpasses
    if the track date is older than the given latest

    :param overpasses: list of track date and track id pairs
    :type overpasses: list
    :param start_date: start date in the format of YYYY-MM-DD
    :type start_date: str
    :return: overpasses that are younger or equal to the start date
    :rtype: list
    """
    latest_date = datetime.strptime(latest_date, "%Y-%m-%d").date()
    return [
        overpass
        for overpass in overpasses
        if datetime.strptime(overpass[0], "%Y-%m-%d").date() >= latest_date
    ]


def get_icesat2_data(
    aoi: gpd.GeoDataFrame,
    clip_geom: bool = True,
    icesat2_product: str = "atl03",
    series_included: list = ["Medium", "High"],
    start: date = date.min,
    stop: date = date.max,
) -> pd.DataFrame:
    """
    Get data from the openaltimetry web API for icesat 2 data.

    Args:
        aoi: geodataframe containing the area of interest.
        clip_geom: whether to clip resulting data from the API by the aoi.
        icesat2_project: product of the icesat api to use.
        series_included: which icesat quality to use.
        start: start filter to use when requesting data from the api.
        stop: stop filter to use when requesting data from the api.

    Returns:
        dataframe containing icesat series.
    """
    buffered_aoi = buffer_geometry(aoi, buffer=5000)
    bounding_box = get_bounding_box(buffered_aoi, icesat2_product)

    df_concat = pd.DataFrame()
    if clip_geom:
        clip_feature = get_clip_feature(aoi)

    for bounds in bounding_box:
        track_ids = get_track_ids(bounds)
        overpasses = get_overpasses(track_ids, start, stop)

        for track_date, track_id in overpasses:
            # This function will request the 6 beams data using OpenAltimetry's API
            photon_data = get_icesat2_data_from_OA_api(
                icesat2_product=icesat2_product,
                boundingbox=bounds,
                track_date=track_date,
                track_id=track_id,
                series_included=series_included,
            )

            if icesat2_product == "atl03":
                rows = parse_atl03(photon_data, track_id, track_date)

            elif icesat2_product == "atl08":
                rows = parse_atl08(photon_data, track_id, track_date)

            if rows:
                logger.info(
                    f"downloaded data for date {track_date} and track id {track_id}"
                )

                df = pd.DataFrame(rows)
                if clip_geom:
                    df_filtered = filter_icesat2_data(df, clip_feature)
                    df = df_filtered

                if df.empty:
                    logger.info(
                        f"no points within reservoir convex hull for track date {track_date} and track id {track_id}"
                    )
                    continue

                df = process_icesat2_for_gee(df, icesat2_product)

                df_concat = pd.concat([df_concat, df])

    if df_concat.empty:
        logger.warning("No icesat2 data found for given aoi")
    return df_concat


def get_latest_overpass_date(geom: gpd.GeoDataFrame) -> date:
    """Retrieves the latest date of an icesat2 overpass for a given geometry

    :param geom: Geometry of area of interest
    :type geom: gpd.GeoDataFrame
    :return: Date of latest icesat2 overpass
    :rtype: datetime.date
    """
    buffered_geom = buffer_geometry(geom, buffer=5000)
    bbox = get_bounding_box(buffered_geom)
    track_ids = get_track_ids(bbox[0])
    overpasses = get_overpasses(track_ids)
    if not overpasses:  # Overpasses can be empty if reservoir is too small
        return None
    latest_date = sorted(
        [datetime.strptime(date_pair[0], "%Y-%m-%d").date() for date_pair in overpasses]
    )[-1]
    return latest_date
