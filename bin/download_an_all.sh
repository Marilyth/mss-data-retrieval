#!/bin/bash
mars <<EOF
retrieve,
    class = od,
    time = $2,
    date = $1,
    levelist=1/to/137,
    levtype=ml,
    param=130.128/131.128/132.128/133.128/135.128/152.128/155.128/203.128/246.128/247.128/248.128,
    stream=oper,
    type=fc,
    area=0/0/-80/360,
    grid=1.0/1.0,
    target="grib/$1T$2.fc.ml.grib"

retrieve,
    levelist=all,
    levtype=pv,
    param=3/54/131/132/133/203,
    target="grib/$1T$2.fc.pv.grib"

retrieve, 
    levelist=all,
    levtype=pt,
    param=54/60/131/132/133/155/203,
    target="grib/$1T$2.fc.tl.grib"

retrieve,
    levtype=sfc,
    param=151.128/165.128/166.128/186.128/187.128/188.128,
    target="grib/$1T$2.fc.sfc.grib"
EOF
