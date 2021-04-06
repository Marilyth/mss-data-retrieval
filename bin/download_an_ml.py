import sys
from ecmwfapi import ECMWFService, ECMWFDataServer

# To run this example, you need an API key
# available from https://api.ecmwf.int/v1/key/

#server = ECMWFService("era15")
server = ECMWFDataServer()

"""server.execute({
    "class": "od",
    "date": sys.argv[1],
    "time": sys.argv[2],
    "expver": 1,
    "levelist": sorted(range(1, 138)),
    "levtype": "ml",
    "param": ["129", "130", "131", "132", "133", "135", "152", "155", "203", "246", "247", "248"],
    "stream": "oper",
    "area": "0/0/-80/360",
    "grid": "1.0/1.0",
    "type": "an"
    },
    "grib/" + sys.argv[1] + "T" + sys.argv[2] + ".an.ml.grib")"""

server.retrieve({
    'dataset' : 'yotc',
    'time'    : sys.argv[2],
    'date'    : sys.argv[1],
    'type'    : "an",
    'levelist': "all",
    'levtype' : "ml",   
    'class'   : "od",
    'area'    : "0/0/-80/360",
    'param'   : ["129", "130", "131", "132", "133", "135", "152", "155", "203", "246", "247", "248"],
    'grid'    : "1.0/1.0",
    'target'  : f"grib/{sys.argv[1]}T{sys.argv[2]}.an.ml.grib"
    })
