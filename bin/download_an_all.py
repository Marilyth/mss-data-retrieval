import sys
from threading import Thread
# To run this example, you need an API key
# available from https://api.ecmwf.int/v1/key/

import cdsapi
c_ml = cdsapi.Client()
c_pl = cdsapi.Client()
c_pv = cdsapi.Client()
c_tl = cdsapi.Client()
c_sfc = cdsapi.Client()

def ml():
    c_ml.retrieve('reanalysis-era5-complete', {
        'class': 'ea',
        'time': sys.argv[2],
        'date': sys.argv[1],
        'expver': '1',
        'levelist': 'all',
        'levtype': 'ml',
        'param': "all",
        'stream': 'oper',
        'type': 'an',
        "area": "0/0/-80/360",
        "grid": "1.0/1.0",
    }, f'grib/era5{sys.argv[1]}T{sys.argv[2]}.ml.grib')

def pl():
    c_pl.retrieve('reanalysis-era5-complete', {
        'class': 'ea',
        'time': sys.argv[2],
        'date': sys.argv[1],
        'expver': '1',
        'levelist': 'all',
        'levtype': 'pl',
        'param': 'all',
        'stream': 'oper',
        'type': 'an',
        "area": "0/0/-80/360",
        "grid": "1.0/1.0",
    }, f'grib/era5{sys.argv[1]}T{sys.argv[2]}.pl.grib')

def pv():
    c_pv.retrieve('reanalysis-era5-complete', {
        'class': 'ea',
        'time': sys.argv[2],
        'date': sys.argv[1],
        'expver': '1',
        'levelist': 'all',
        'levtype': 'pv',
        'param': 'all',
        'stream': 'oper',
        'type': 'an',
        "area": "0/0/-80/360",
        "grid": "1.0/1.0",
    }, f'grib/era5{sys.argv[1]}T{sys.argv[2]}.pv.grib')

def tl():
    c_tl.retrieve('reanalysis-era5-complete', {
        'class': 'ea',
        'time': sys.argv[2],
        'date': sys.argv[1],
        'expver': '1',
        'levelist': [340, 360],
        'levtype': 'pt',
        'param': 'all',
        'stream': 'oper',
        'type': 'an',
        "area": "0/0/-80/360",
        "grid": "1.0/1.0",
    }, f'grib/era5{sys.argv[1]}T{sys.argv[2]}.tl.grib')

def sfc():
    c_sfc.retrieve('reanalysis-era5-complete', {
        'class': 'ea',
        'time': sys.argv[2],
        'date': sys.argv[1],
        'expver': '1',
        'levelist': 'all',
        'levtype': 'sfc',
        'param': 'all',
        'stream': 'oper',
        'type': 'an',
        "area": "0/0/-80/360",
        "grid": "1.0/1.0",
    }, f'grib/era5{sys.argv[1]}T{sys.argv[2]}.sfc.grib')

threads = []
#threads.append(Thread(target=ml))
#threads.append(Thread(target=pl))
#threads.append(Thread(target=pv))
threads.append(Thread(target=tl))
#threads.append(Thread(target=sfc))

for thread in threads:
    thread.start()

for thread in threads:
    thread.join()
