import io
import os
import re

from setuptools import find_packages
from setuptools import setup


def read(filename):
    filename = os.path.join(os.path.dirname(__file__), filename)
    text_type = type(u"")
    with io.open(filename, mode="r", encoding='utf-8') as fd:
        return re.sub(text_type(r':[a-z]+:`~?(.*?)`'), text_type(r'``\1``'), fd.read())


setup(
    name="bathypy",
    version="0.1.0",
    description="BAthymetry estimation Through HYdrological Remote Sensing with Python",
    long_description=read("README.rst"),
    url="https://github.com/global-water-watch/bathypy",
    license='MIT',
    author="Hessel Winsemius",
    author_email="hessel.winsemius@deltares.nl",
    packages=find_packages(exclude=('tests',)),
    package_dir={"bathypy": "bathypy"},
    python_requires=">=3.7",
    install_requires=[
        "aiohttp",
        "cartopy",
        "descartes",
        "earthengine-api",
        "eepackages",
        "fastparquet",
        "gcsfs",
        "geemap",
        "geojson",
        "geopandas",
        "google-cloud-storage",
        "hydromt",
        "jupyter",
        "matplotlib",
        "netCDF4",
        "numba",
        "numpy",
        "pip",
        "pyarrow",
        "pyproj",
        "pytest",
        "rasterio",
        "requests",
        "retry",
        "scipy",
        "tqdm",
        "xarray",
    ],
    include_package_data=True,
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
)
