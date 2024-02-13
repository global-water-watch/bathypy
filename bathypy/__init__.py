__version__ = "0.1.0"

# ensure to use pygeos
import os
os.environ['USE_PYGEOS'] = '0'
from bathypy import storage
from bathypy import icesat2
from bathypy import gww_requests
from bathypy import static_hypso
