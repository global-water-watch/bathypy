from logging.config import dictConfig
import os
from pathlib import Path
import random

import ee
import geopandas as gpd
import pytest
from shapely.geometry import shape
import xarray as xr
from yaml import safe_load

from bathypy.static_hypso import Bathymetry
# from bathypy.surface import get_waterbody
from bathypy.utils import gee_feature_to_gdf
from bathypy.icesat2 import get_bounding_box
from bathypy import gww_requests

EXAMPLE_DATA_DIR = os.path.join(os.path.split(__file__)[0], "..", "examples")


@pytest.fixture(scope="session")
def ee_init():
    ee.Initialize(opt_url="https://earthengine-highvolume.googleapis.com")
    os.environ["GCLOUD_PROJECT"] = "global-water-watch"


def get_reservoir_id() -> int:
    with open(Path(__file__).absolute().parent / "data/test_reservoirs.txt", "r") as f:
        reservoir_ids = f.readlines()
        return int(random.choice(reservoir_ids).strip("\n"))


@pytest.fixture
def fn_topography():
    return os.path.join(EXAMPLE_DATA_DIR, "datasets", "topo_mita.nc")


@pytest.fixture
def polygon():
    reservoir_id = 88643  # Mita Hills, Zambia. Or any other ID.
    response = gww_requests.get_reservoir(reservoir_id)
    return shape(response.json()["geometry"])


@pytest.fixture
def ds_topography(fn_topography):
    return xr.open_dataset(fn_topography)


@pytest.fixture
def bathymetry(polygon):
    return Bathymetry(polygon=polygon)


@pytest.fixture
def bathymetry_with_data(bathymetry, ds_topography):
    bathymetry.set_topography(ds_topography)
    return bathymetry


@pytest.fixture
def bathymetry_interpolated(bathymetry_with_data):
    bathymetry_with_data.skeletonize_bathymetry(iter=2)
    bathymetry_with_data.interpolate()
    return bathymetry_with_data


@pytest.fixture(scope="session")
def uses_logging():
    with open(Path(__file__).absolute().parent.parent / "logging.yaml", "r") as f:
        config = safe_load(f.read())
        dictConfig(config)


# @pytest.fixture(scope="session")
# def feature() -> ee.Feature:
#     fid = get_reservoir_id()
#     return get_waterbody(fid)


# @pytest.fixture(scope="session")
# def aoi() -> gpd.GeoDataFrame:
#     fid = get_reservoir_id()
#     return gee_feature_to_gdf(get_waterbody(fid), fid)


@pytest.fixture(scope="session")
def bbox(aoi: gpd.GeoDataFrame) -> list:
    return get_bounding_box(aoi, icesat2_product="atl03")