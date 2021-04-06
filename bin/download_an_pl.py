import sys
from ecmwfapi import ECMWFDataServer

# To run this example, you need an API key
# available from https://api.ecmwf.int/v1/key/

server = ECMWFDataServer()

server.retrieve({
    "class": "yt",
    "dataset": "yotc",
    "date": "2010-04-01",
    "expver": "1",
    "levelist": "1/2/3/5/7/10/20/30/50/70/100/150/200/250/300/400/500/600/700/800/850/900/925/950/1000",
    "levtype": "pl",
    "param": "60.128/129.128/130.128/131.128/132.128/133.128/135.128/138.128/155.128/157.128/203.128",
    "step": "0",
    "stream": "oper",
    "time": "12:00:00",
    "type": "an",
    "target": "grib/test.an.pl.grib",
})

"""server.retrieve({
    'dataset' : 'yotc',
    'time'    : sys.argv[2],
    'date'    : sys.argv[1],
    'type'    : "an",
    'levelist': "1/2/3/5/7/10/20/30/50/70/100/150/200/250/300/400/500/600/700/800/850/900/925/950/1000",
    'levtype' : "pl",   
    'class'   : "od",
    'area'    : "0/0/-80/360",
    'param'   : ["129", "130", "131", "132", "133", "135", "152", "155", "203", "246", "247", "248"],
    'grid'    : "1.0/1.0",
    'target'  : f"grib/{sys.argv[1]}T{sys.argv[2]}.an.pl.grib"
    })"""
