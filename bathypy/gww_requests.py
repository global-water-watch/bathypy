import geopandas as gpd
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import cartopy
import cartopy.io.img_tiles as cimgt
import cartopy.crs as ccrs
import datetime

from dateutil.relativedelta import relativedelta
from scipy import stats
from shapely.geometry import shape, box, mapping

base_url = "https://api.globalwaterwatch.earth"


def to_geopandas(data):
    """
    Ingests list of reservoirs and converts into a geopandas GeoDataFrame for further analyses

    """
    geoms = [shape(f["geometry"]) for f in data]
    props = [{**f["properties"], **{"id": f["id"]}} for f in data]
    return gpd.GeoDataFrame(props, geometry=geoms, crs=4326)


def get_reservoirs(skip=1, limit=5, base_url=base_url):
    """
    Gets reservoirs from API. Return dict with IDs.

    """
    url = f"{base_url}/reservoir"
    params = {
        "skip": skip,
        "limit": limit,
    }
    return requests.get(url, params=params)


def get_reservoir(reservoir_id):
    """
    Get reservoir (geometry and props) by ID
    """
    url = f"{base_url}/reservoir/{reservoir_id}"
    return requests.get(url)


def get_reservoirs_by_geom(geom, base_url=base_url):
    """
    Gets reservoirs from API. Return dict with IDs.

    """
    url = f"{base_url}/reservoir/geometry"
    # do post request to end point with the serialized geometry as post data
    return requests.post(url, data=geom)


def get_reservoir_ts(reservoir_id, start, stop):
    """
    Get time series data for reservoir with given ID
    """
    url = f"{base_url}/reservoir/{reservoir_id}/ts"
    params = {
        "start": start.strftime("%Y-%m-%dT%H:%M:%S"),
        "stop": stop.strftime("%Y-%m-%dT%H:%M:%S")
    }
    return requests.get(url, params=params)


def plot_features_map(feats, ax=None, figsize=(20, 13), tiles=None, zoom_level=1, tiles_kwargs={}, **kwargs):
    """
    add a set of features to a GeoAxes map
    """
    if ax is None:
        f = plt.figure(figsize=figsize)
        if tiles is not None:
            tiler = getattr(cimgt, tiles)(**tiles_kwargs)
            crs = tiler.crs
        else:
            crs = ccrs.PlateCarree()
        # make point collection
        ax = plt.subplot(projection=crs)
        if tiles is not None:
            ax.add_image(tiler, zoom_level, zorder=1)
            feats.to_crs(3857).plot(ax=ax, zorder=2, **kwargs)
        else:
            feats.plot(ax=ax, **kwargs)
    return ax