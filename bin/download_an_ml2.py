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
    'levelist': 'all',
    'levtype': 'ml',
    'param': ["129", "130", "131", "132", "133", "135", "152", "155", "203", "246", "247", "248"],
    'stream': 'oper',
    'type': 'an',
}, f'grib/era5{sys.argv[1]}T{sys.argv[2]}.ml.grib')

