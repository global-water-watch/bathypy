Surface2Storage - Global Water Watch 
====================================

.. image:: https://img.shields.io/pypi/v/gwwstorage.svg
    :target: https://pypi.python.org/pypi/gwwstorage
    :alt: Latest PyPI version

This package provides several approaches for derivation of estimated hypsometry (i.e. relationships between water level, surface area and stored volume in lakes and reservoirs) from earth observations. Two general methods are available:
1) derivation of hypsometry from a combination of surface area estimates and ICESat-2 water level estimates, paired on dates.
2) derivation of geostatistically inference of bathyemtry using static terrain datasets and river topology.

This repository has been made possible thanks to grant 4000139266 of the European Space Agency's EO4Society program. Deltares and the NEtherlands Red Cross are very grateful for this contribution.

Installation
============

Prerequisites
-------------

BathyPy makes use of a number of cloud platforms to perform the processing. In detail the following is needed to run the methods in this repository.

* A Google Earth Engine (GEE) account, needed to perform large-scale satellite data analysis for surface water estimation.
* A Google Cloud Storage: this is needed to store results from the GEE analyses and store intermediate ICESat-2 retrievals, retrieved from the ICESat-2 API, enabling its use in GEE.

Setting up a Google Earth Engine account
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To set up a google earth engine account, please follow the instructions on https://earthengine.google.com/new_signup/

Setting up a Google Cloud Storage Bucket
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To set up a Google Cloud Storage, please follow the instructions on https://cloud.google.com/storage. Here click on "Go to console" to make a new google cloud bucket.

Installation of this package
----------------------------

To run the scripts and notebooks in this repository, you can follow the instructions below to get a clean conda-based
virtual environment.

First, clone the code with `git` and move into the cloned folder.

.. code-block:: console

    git clone git@github.com:global-water-watch/bathypy.git
    cd bathypy

Setup a virtual environment as follows:

.. code-block:: console

    conda env create -f environment.yml

Connecting python api to cloud storage
--------------------------------------
To make use of the cloud storage of Global Water Watch, you will need a private key.
Follow these instructions to obtain one and set it. Your Google account must have been added to the Water Watch project
before you can access this storage.

Go to https://console.cloud.google.com/ and login if necessary.
Ensure at the very top of the page, that your cloud storage project is selected.
Below we assume your cloud storage bucket is called `Water Watch`.

In the left menu, go to IAM & Admin --> Service Accounts. In the table, click on
`storage@water-watch.iam.gserviceaccount.com` or similar.

In the appearing page, select `KEYS`, and then `ADD KEY` --> `Create new key`.
Select JSON as key type and then select `CREATE`. A JSON file will be downloaded.

Now store this JSON file in a easy to locate folder. You can place it under this repos in a folder
`keys`. We will assume this in the rest of this instruction.

Now add an environment variable to your system environment variables as follows:

.. code-block:: console

    GOOGLE_APPLICATION_CREDENTIALS=/path/to/keyfile.json

Now you can use your storage bucket, the easiest approach to use it is to use the function
`gwwstorage.cloud.get_bucket`.

After this is set, you can install the package. To install it as a user, please type.

.. code-block:: console

    pip install .

To install the package for development, clone and run:

.. code-block:: console

    pip install -e .

Notebooks
---------

To run the notebooks, use `conda` together with the `environment.yml`.

Docker
------
You can also run in Docker. Note that it has to build the conda environment from scratch, so this can take up to
10 minutes to complete.

.. code-block:: console

    docker build . -f docker/notebooks/Dockerfile -t bathypy-nb
    docker run -p 8888:8888 -v ~/.config:/home/jovyan/.config -v $(pwd):/home/jovyan/work bathypy-nb

Background and examples of methods
==================================

Please refer to our documentation pages for guidance and backgrounds of the methods. In this
documentation we have included 3 example notebooks that provide the main features.

Background of method
~~~~~~~~~~~~~~~~~~~~

Lorem ipsum


Retrieving and uploading ICESat-2 tracks over polygon
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Icesat2 data can be retrieved for a polygon in the form of a geopandas.GeoDataFrame
in EPSG 4326. The data is returned as a pandas.DataFrame with point coordinates as
lon, lat columns.

.. code-block:: python

    import geopandas as gpd
    from gwwstorage.altimetry import get_dynamic_elevation_data
    
    aoi = gpd.GeoDataFrame()
    icesat2_data = get_dynamic_elevation_data(
        aoi=aoi,
        dataset="icesat2",
        icesat2_product="atl03"
        )


    ...

The icesat2 data can be uploaded to a google cloud bucket

.. code-block:: python

    from gwwstorage.altimetry import cloud_store_dynamic_elevation_data

    cloud_store_dynamic_elevation_data(
        df=icesat2_data,
        blob_name="filename",
        bucket_name="global-water-watch-cache",
        gcs_folder="icesat2-cache"
        )
    ...

Establishing surface area estimates at ICESat-2 overpass dates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lorem ipsum

Retrieving paired surface area and water level series
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lorem ipsum


Bathymetry and hypsometry points from static elevation datasets
---------------------------------------------------------------

Lorem ipsum


Requirements
------------

Compatibility
-------------

Licence
-------

