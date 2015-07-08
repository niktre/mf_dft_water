#!/bin/bash
eps=$1
Lx=$2
Ly="1.6"
Lz="16.4"
Cx=$3
Cy=${Ly}
Cz="2."
den="33.40"
lambda=$4
restart=$5
nDims="2"
NVT="0"
qsub -v den=${den},eps=${eps},lambda=${lambda},restart=${restart},Lx=${Lx},Ly=${Ly},Lz=${Lz},Cx=${Cx},Cy=${Cy},Cz=${Cz},nDims=${nDims},NVT=${NVT} submit_gc_groove.sh
