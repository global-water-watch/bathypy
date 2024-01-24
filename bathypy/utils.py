from typing import Sequence

import numba
from numba import jit, njit
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
from shapely import unary_union
import ee
import geopandas as gpd


@njit(parallel=True)
def points_in_polygon_parallel(
    coordinates: np.ndarray, polygon: np.ndarray
) -> np.ndarray:
    """This function uses a parallelized approach for determining whether multiple
    points are within a polygon.

    :param coordinates: list of coordinate pairs of points
    :type coordinates: np.ndarray
    :param polygon: collection of points forming a polygon
    :type polygon: np.ndarray
    :return: A boolean array with the same length of the list point coordinates
    :rtype: np.ndarry
    """

    ln = len(coordinates)
    D = np.empty(ln, dtype=numba.boolean)
    for i in numba.prange(ln):
        coords = coordinates[i]
        D[i] = point_in_polygon(polygon, coords)
    return D


@jit(nopython=True)
def point_in_polygon(polygon: np.ndarray, point: np.ndarray) -> bool:
    length = len(polygon)
    intersections = 0

    dx2 = point[0] - polygon[0][0]
    dy2 = point[1] - polygon[0][1]
    jj = 1

    while jj < length:
        dx = dx2
        dy = dy2
        dx2 = point[0] - polygon[jj][0]
        dy2 = point[1] - polygon[jj][1]

        F = (dx - dx2) * dy - dx * (dy - dy2)
        if 0.0 == F and dx * dx2 <= 0 and dy * dy2 <= 0:
            return True

        if (dy >= 0 and dy2 < 0) or (dy2 >= 0 and dy < 0):
            if F > 0:
                intersections += 1
            elif F < 0:
                intersections -= 1

        jj += 1

    return intersections != 0


def get_bbox_from_feature(feature: ee.Feature) -> tuple:
    """Retrieves a bounding box from a earth engine feature

    :param feature: feature containing a geometry with projection
    :type feature: ee.feature.Feature
    :return: a tuple of min and max xy coordinates
    :rtype: tuple
    """
    feature_bounds = feature.geometry().bounds().getInfo()

    geom = feature_bounds["coordinates"][0]
    minx, miny, maxx, maxy = Polygon(geom).bounds
    return (minx, miny, maxx, maxy)


def bounding_box_tiles(bbox_coords: Sequence[float], step_size: float = 0.5) -> list:
    """Divides a bounding box to smaller bounding boxes given a stepsize

    :param bbox_coords: The min x, min y, max x, and max y of a bounding box
    :type bbox_coords: tuple or list
    :param step_size: The length of the smaller bounding boxes
    :type step_size: float
    :return: list of bounding boxes with min x, min y, max x, max y coordinates
    :rtype: list
    """
    minx, miny, maxx, maxy = bbox_coords
    box_list = []
    for y in np.arange(miny, maxy, step_size):
        for x in np.arange(minx, maxx, step_size):
            xmin = x
            ymin = y
            xmax = x + step_size
            ymax = y + step_size
            box = [xmin, ymin, xmax, ymax]
            box_list.append(box)
    return box_list


def gee_feature_to_gdf(feature: ee.Feature, fid: int) -> gpd.GeoDataFrame:
    """Parses a Google Earth Engine feature to Geopandas GeoDataFrame.

    :param feature: An earth engine feature with a polygon geometry.
    :type feature: ee.Feature
    :param fid: Id of the feature
    :type fid: int
    :return:  A dataframe containing the geometry and fid of the feature.
    :rtype: gpd.GeoDataFrame
    """
    # print(fid)
    geoms = feature.getInfo()["geometry"]["coordinates"]

    ### Note Antonio: Some reservoirs are Multiple polys, This is an attempt to merge them prior to the clipping of Icesat2 points (unary_union).
    # If elements are disconnected this will fail, in such case I revert to get the largest polygon.

    if len(geoms) > 1:
        try:
            polygon = unary_union(
                MultiPolygon([Polygon(coords) for coords in geoms])
            )

        except:
            polygon = unary_union(
                MultiPolygon([Polygon(coords[0]) for coords in geoms])
            )

    else:
        polygon = Polygon(geoms[0])

    if isinstance(polygon, MultiPolygon): ## If it is still a polygon, because the unary union does not resolve a joint geometry, get the largest.
            polygon = max((poly for poly in polygon.geoms), key=lambda poly: poly.area)


    gdf = gpd.GeoDataFrame(data={"fid": [id], "geometry": [polygon]}, crs=4326)

    return gdf


def buffer_geometry(aoi: gpd.GeoDataFrame, buffer: int) -> gpd.GeoDataFrame:
    """Buffer the first geometry of a geopandas GeoDataFrame.

    :param aoi: A GeoDataFrame containing a geometry
    :type aoi: gpd.GeoDataFrame
    :param buffer: buffer distance in meters
    :type buffer: int
    :return: A GeoDataFrame with a buffered geometry
    :rtype: gpd.GeoDataFrame
    """
    geoseries = aoi.loc[:, "geometry"]
    proj_crs = geoseries.estimate_utm_crs()
    gdf = aoi.to_crs(proj_crs)
    gdf.loc[0, "geometry"] = gdf.loc[0, "geometry"].buffer(buffer)
    gdf = gdf.to_crs(4326)
    return gdf
