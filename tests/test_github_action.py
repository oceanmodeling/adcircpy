from pydap.client import open_url
from netCDF4 import Dataset
import xarray as xr

url = 'https://icdc.cen.uni-hamburg.de/thredds/dodsC/ftpthredds/hamtide//m2.hamtide11a.nc'

test_get_pydap():
    ds = open_url(url)
    print("OPENDAP", ds)

test_open_w_nc4():
    ds = Dataset(url)
    print("NETCDF4", ds)

test_open_w_xr():
    ds = xr.open_dataset(url)
    print("XARRAY", ds)

