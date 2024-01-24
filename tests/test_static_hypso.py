import numpy as np

from gwwstorage.static_hypso import Bathymetry
from shapely.geometry import Polygon, MultiPolygon
from xarray import Dataset


def test_polygon(polygon):
    assert(isinstance(polygon, Polygon) or isinstance(polygon, MultiPolygon))

def test_fn(fn_topography):
    print(fn_topography)

def test_topo(ds_topography):
    # ensure we are looking at a dataset with topo
    assert(isinstance(ds_topography, Dataset))


def test_bathymetry(bathymetry, ds_topography):
    # b = bathymetry_with_data
    bathymetry.set_topography(ds_topography)
    assert(isinstance(bathymetry.polygon, Polygon) or isinstance(bathymetry.polygon, MultiPolygon))


def test_bathymetry_with_data(bathymetry_with_data):
    # ensure we are looking at a dataset with topo
    assert(isinstance(bathymetry_with_data, Bathymetry))


def test_skeletonize(bathymetry_with_data):
    bathymetry_with_data.skeletonize_bathymetry(iter=2)

def test_interpolate(bathymetry_with_data):
    bathymetry_with_data.skeletonize_bathymetry(iter=2)
    bathymetry_with_data.interpolate()
    # resulting bathymetry should not contain any NaN
    assert(not(np.isnan(bathymetry_with_data.ds_topography["bathymetry"].values).any()))

def test_hypso(bathymetry_interpolated):
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use("Qt5Agg")
    hypsometry = bathymetry_interpolated.get_hypsometry(add_to_max_h=10)
    bathymetry_interpolated.plot_topo(blend_mode="soft", vert_exag=40)
    # ax = hypsometry.plot()
    plt.show()
