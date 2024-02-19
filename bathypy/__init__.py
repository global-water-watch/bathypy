__version__ = "0.1.0"

# ensure to use pygeos
import os
os.environ['USE_PYGEOS'] = '0'

# older versions of rasterio with ogr have issues importing versions, here, fix the required imports
try:
    import osgeo
    from rasterio._version import gdal_version, get_geos_version, get_proj_version
except:
    # osgeo no longer needed with newer versions of rasterio
    pass

from bathypy import storage
from bathypy import icesat2
from bathypy import gww_requests
from bathypy import static_hypso
