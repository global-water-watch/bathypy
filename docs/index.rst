.. BathyPy documentation master file, created by
   sphinx-quickstart on Fri Nov 19 15:06:11 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _index_page:

==========================================
Welcome to BathyPy's documentation!
==========================================

**BathyPy**, is a fully Open Source library for estimating bathymetry from both static and dynamic space-borne
datasets of topography and water surface areas. It builds upon Earth Observation methods of
`Global Water Watch <https://globalwaterwatch.earth>`_ a Google funded project that measures and discloses water storage
in reservoirs and lakes world wide.

This Open source library encapsulates methods to derive the bathymetry and hypsometric relationships of lakes and
reservoirs and consists of two different methods, both starting with a broadly defined polygon that encapsulates
the lake or reservoir surface, and providing a set of commensurate water levels and surface areas. These can then be used
to construct a parametric relationship between water levels and surface areas or (when integrated over water level)
volumes. The methods are in a nutshell:

+----------------------------------+-----------------------------------------------------------------------------------+-----------------+
| Method                           | description                                                                       | Applicable size |
+==================================+===================================================================================+=================+
| Dynamic bathymetry               | Uses a time series of dynamically estimated surface areas in a polygon of interest| >= 3 km3        |
| estimation                       | and time series of ATL03 or ATL08 level processed altimetry from the ICESat-2     |                 |
|                                  | altimetry mission. It pairs the water levels and surface area using the closest   |                 |
|                                  | matching date between the two.                                                    |                 |
+----------------------------------+-----------------------------------------------------------------------------------+-----------------+
| Static bathymetry                | uses gridded data from a static topography mission such as SRTM                   | < 3 km3         |
| estimation                       | along with derived products, including the Digital Elevation Model (DEM)          |                 |
|                                  | D8 flow directions, upstream surface area, strahker stream order and slope along  |                 |
|                                  | the stream to detect areas in the DEM that likely are water surface, maks these.  |                 |
|                                  | Uses geomorphological relationships to interpolate bathymetry following the       |                 |
|                                  | dominant flow paths.                                                              |                 |
+----------------------------------+-----------------------------------------------------------------------------------+-----------------+

The dynamic methods rely heavily on Google Earth Engine and associated cloud storage
and therefore require a somewhat demanding installation procedure to prepare. The static bathymetry methods rely on the
`HydroMT <https://deltares.github.io/hydromt>`_ package.

For whom?
=========
**BathyPy** provides essential information on lakes and reservoirs and has many potential use cases including:

- Water and dam managers, that wish to estimate the water volume in reservoirs from local water level observations or
  space-based water level or water surface areas.
- Hydrologists, and hydrological modelling for water resource planning and allocation.
- Basin planners in basins with shared water resources with upstream countries, who have little or no information on
  upstream stored water.
- Scientists in the field of hydrology, ecology, lymnology and climate.

.. note::

   Acknowledgement: This software has been created under a ESA grant with contract no. 4000139266/22/I-DT-bgh, Future EO-1
   EO Science for Society, project name: "Surface-2-Storage: Enhancing drought Resilience by global real-time monitoring
   of water volumes stored in small to medium sized reservoirs."

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Table of Contents

   installation
   quickstart
   api.rst

