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
        'levelist': sorted(range(1, 138)),
        'levtype': 'ml',
        'param': [129.128, 130.128, 131.128, 132.128, 133.128, 135.128, 152.128, 155.128, 203.128, 246.128, 247.128, 248.128],
        'stream': 'oper',
        'type': 'an',
        "area": "90/0/-90/360",
        "grid": "1.0/1.0",
    }, f'grib/{sys.argv[1]}T{sys.argv[2]}.an.ml.grib')

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
        "area": "90/0/-90/360",
        "grid": "1.0/1.0",
    }, f'grib/{sys.argv[1]}T{sys.argv[2]}.an.pl.grib')

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
        "area": "90/0/-90/360",
        "grid": "1.0/1.0",
    }, f'grib/{sys.argv[1]}T{sys.argv[2]}.an.pv.grib')

def tl():
    c_tl.retrieve('reanalysis-era5-complete', {
        'class': 'ea',
        'time': sys.argv[2],
        'date': sys.argv[1],
        'expver': '1',
        'levelist': 'all',
        'levtype': 'pt',
        'param': 'all',
        'stream': 'oper',
        'type': 'an',
        "area": "90/0/-90/360",
        "grid": "1.0/1.0",
    }, f'grib/{sys.argv[1]}T{sys.argv[2]}.an.tl.grib')

def sfc():
    c_sfc.retrieve('reanalysis-era5-complete', {
        'class': 'ea',
        'time': sys.argv[2],
        'date': sys.argv[1],
        'expver': '1',
        'levtype': 'sfc',
        'param': ["151.128", "165.128", "166.128", "186.128", "187.128", "188.128"],
        'stream': 'oper',
        'type': 'an',
        "area": "90/0/-90/360",
        "grid": "1.0/1.0",
    }, f'grib/{sys.argv[1]}T{sys.argv[2]}.an.sfc.grib')

threads = []
if "ml" in sys.argv:
    threads.append(Thread(target=ml))
if "pl" in sys.argv:
    threads.append(Thread(target=pl))
if "pv" in sys.argv:
    threads.append(Thread(target=pv))
if "tl" in sys.argv:
    threads.append(Thread(target=tl))
if "sfc" in sys.argv:
    threads.append(Thread(target=sfc))

for thread in threads:
    thread.start()

for thread in threads:
    thread.join()
