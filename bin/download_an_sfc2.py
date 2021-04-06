import sys

# To run this example, you need an API key
# available from https://api.ecmwf.int/v1/key/

import cdsapi
c = cdsapi.Client()
c.retrieve('reanalysis-era5-complete', {
    'class': 'ea',
    'time': sys.argv[2],
    'date': sys.argv[1],
    'expver': '1',
    'levtype': 'sfc',
    'param': ["151", "165", "166", "186", "187", "188"],
    'stream': 'oper',
    'area': "0/0/-80/360",
    'type': 'an',
    'grid': "1.0/1.0",
}, f"grib/{sys.argv[1]}T{sys.argv[2]}.an.sfc.grib")

