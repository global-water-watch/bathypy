.. _installation:

============
Installation
============

.. _user_install:

User install
============

To install **BathyPy**, you will need a package manager in the Anaconda/Miniconda ecosystem such as **conda** or
**mamba**.

We recommend using the Mambaforge_ Python distribution. This installs Python and the mamba package manager. 
Miniforge_ and Miniconda_ will install Python and the conda package manager.

In general, **mamba** is a lot faster than **conda** and requires less memory.

.. _install_conda-forge:

Get code from github
====================
Retrieve the latest code from github as follows:

.. code-block:: console

    $ git clone https://github.com/global-water-watch/research-reservoir-storage
    $ cd research-reservoir-storage

Setup of conda environment
==========================
We always recommend working in a dedicated virtual environment for **BathyPy**. We have a fully prepared environment
file that installs all required dependencies as follows:

.. code-block:: console

    $ conda env create -f environment.yml
    $ conda activate BathyPy

Once the environment is prepared, install the package as follows:

.. code-block:: console

    $ pip install .

If you wish to develop in the code base then install as follows:

.. code-block:: console

    $ pip install -e .

Connecting Google Earth Engine
------------------------------

If you wish to use the dynamic methods, then Google Earth Engine and a Google Cloud Bucket are needed.

Setting up a Google Earth Engine account
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To set up a google earth engine account, please follow the instructions on https://earthengine.google.com/new_signup/

Setting up a Google Cloud Storage Bucket
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To set up a Google Cloud Storage, please follow the instructions on https://cloud.google.com/storage. Here click on "Go to console" to make a new google cloud bucket.

Connecting python api to cloud storage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To make use of your cloud storage, you will need a private key to access that cloud storage within BathyPy.
Follow these instructions to obtain one and set it.

Go to https://console.cloud.google.com/ and login if necessary.
Ensure at the very top of the page, that your cloud storage project is selected.
Below we assume your cloud storage bucket is called `Water Watch`.

In the left menu, go to IAM & Admin --> Service Accounts. In the table, click on
`storage@water-watch.iam.gserviceaccount.com` or similar.

In the appearing page, select `KEYS`, and then `ADD KEY` --> `Create new key`.
Select JSON as key type and then select `CREATE`. A JSON file will be downloaded.

Now store this JSON file in an easy to locate folder. You can place it under this repos in a folder
`keys`.

Now add an environment variable to your system environment variables as follows:

.. code-block:: console

    GOOGLE_APPLICATION_CREDENTIALS=/path/to/keyfile.json

Now you can use your storage bucket, the easiest approach to test it is to use the function
`gwwstorage.cloud.get_bucket`.


.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _Mambaforge: https://github.com/conda-forge/miniforge#mambaforge
.. _Miniforge: https://github.com/conda-forge/miniforge
.. _mamba package manager: https://github.com/mamba-org/mamba
.. _conda package manager: https://docs.conda.io/en/latest/
.. _pip package manager: https://pypi.org/project/pip/
.. _manage environments: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
