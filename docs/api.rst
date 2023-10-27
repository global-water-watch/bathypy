.. currentmodule:: gwwstorage.static_hypso

.. _api:

=============
API reference
=============

**BathyPy**'s classes and methods are strongly based on xarray_'s and geopandas_'s data models and
dimensions. The dynamic approach furthermore strongly relies on Google Earth Engine through the ``ee``
python API.

Static methods
==============

.. note::
    For the static methods, we highly recommend consulting the xarray_ user guide side-by-side with this user guide
    to understand how to visualize and interrogate the topography datasets.

In the remaining sections, we describe the API classes, and the functions they are based on.

.. _Bathymetry:

Bathymetry class
================

Class and properties
--------------------

.. autosummary::
    :toctree: _generated

    Bathymetry
    Bathymetry.x
    Bathymetry.y
    Bathymetry.cell_area
    Bathymetry.mask_stream
    Bathymetry.gdf_polygon
    Bathymetry.downstream_basin
    Bathymetry.upstream_area

Setting of properties and attributes
------------------------------------

.. autosummary::
    :toctree: _generated

Adding topography
-----------------

After initializing a ``Bathymetry`` object, a dataset with several topography layers should be set on the ``Bathymetry``
object. Several methods are available.

.. autosummary::
    :toctree: _generated

    Bathymetry.get_topography_from_hydromt
    Bathymetry.read_topography
    Bathymetry.set_topography
    Bathymetry.set_mask_topography
    Bathymetry.set_flow_direction


Interpolating bathymetry
------------------------
The most comprehensive part of ``Bathymetry`` is to fill in missing elevation in lakes. The following methods
provide the necessary functionality

.. autosummary::
    :toctree: _generated

    Bathymetry.skeletonize_bathymetry
    Bathymetry.interpolate

Plotting
--------
Some convenience methods for plotting are available.

.. autosummary::
    :toctree: _generated

    Bathymetry.plot_2d
    Bathymetry.plot_along_stream
    Bathymetry.plot_topo

Hypsometry
----------
After interpolating the bathymetry, a hypsometric relationship in the form of a ``Hypsometry`` object can very easily be
derived.

.. autosummary::
    :toctree: _generated

    Bathymetry.get_hypsometry


Hypsometry class
================
Once a ``Hypsometry`` object is derived, several functionalities are available to fit relationships, work with
the data and export these to figures and files

Class and properties
--------------------

.. autosummary::
    :toctree: _generated

    Hypsometry.water_level
    Hypsometry.area
    Hypsometry.volume
    Hypsometry.area_fit
    Hypsometry.volume_fit
    Hypsometry.power_law_params

Relationships
-------------

.. autosummary::
    :toctree: _generated

    Hypsometry.area_from_wl
    Hypsometry.wl_from_area
    Hypsometry.volume_from_wl
    Hypsometry.volume_from_area

Plotting
--------

.. autosummary::
    :toctree: _generated

    Hypsometry.plot


Dynamic methods
===============

.. note::
    For the dynamic methods, we highly recommend consulting the earth_engine_ user guide side-by-side where it concerns
    computations with Google Earth Engine. Currently, the dynamic methods are applied using several functions.
    We do not yet have a set of classes.

.. currentmodule:: gwwstorage


Surface area detection
----------------------
The functions under subpackage ``surface`` perform water surface area detection in earth_engine_ organize results in
time series form.

.. autosummary::
    :toctree: _generated

    surface.get_surface_area
    surface.get_waterbody
    surface.get_waterbody_ids

Storage
-------

.. autosummary::
    :toctree: _generated

    storage.get_waterbody_images
    storage.compute_surfacewaterarea_singleimage
    storage.compute_surfacewaterarea_JRC
    storage.compute_filled_waterarea
    storage.time_distance
    storage.estimate_watermask_icesatrack
    storage.get_watermasks_icesat2
    storage.sample_from_DEM_datasets
    storage.collect_reservoir_points
    storage.get_watermasks_icesat2_remote


ICESat-2 data retrieval
-----------------------

.. autosummary::
    :toctree: _generated

    icesat2.parse_atl03
    icesat2.parse_atl08
    icesat2.get_bounding_box
    icesat2.get_clip_feature
    icesat2.filter_icesat2_data
    icesat2.process_icesat2_for_gee
    icesat2.drop_redundant_overpasses
    icesat2.get_icesat2_data
    icesat2.get_latest_overpass_date


.. _xarray: https://docs.xarray.dev/
.. _geopandas: https://geopandas.org/
.. _earth_engine: https://developers.google.com/earth-engine/tutorials/community/intro-to-python-api/
