import sys
from ecmwfapi import ECMWFService, ECMWFDataServer

# To run this example, you need an API key
# available from https://api.ecmwf.int/v1/key/

#server = ECMWFService("data")
server = ECMWFDataServer()

"""server.execute({
    "class": "od",
    "date": sys.argv[1],
    "expver": 1,
    "levtype": "sfc",
    "param": ["151", "165", "166", "186", "187", "188"],
    "stream": "oper",
    "area": "0/0/-80/360",
    "grid": "0.25/0.25",
    "time": sys.argv[2],
    "type": "an"
    },
    "grib/" + sys.argv[1] + "T" + sys.argv[2] + ".an.sfc.grib")"""

server.retrieve({
    'dataset' : "yotc",
    'time'    : sys.argv[2],
    'date'    : sys.argv[1],
    'type'    : "an",
    'levtype' : "sfc",
    'stream'  : "oper",
    'param'   : ["151", "165", "166", "186", "187", "188"],#'param'   : "165.128/41.128",
    'area'    : "0/0/-80/360",
    'class'   : "od",
    'grid'    : "1.0/1.0",
    'target'  : f"grib/{sys.argv[1]}T{sys.argv[2]}.an.sfc.grib"
    })
