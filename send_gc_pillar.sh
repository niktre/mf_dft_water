#!/bin/bash
eps=$1
Lx=$2
Ly=${Lx}
Lz="16.4"
Cx=$3
Cy=${Cx}
Cz="2."
den="33.40"
lambda=$4
restart=$5
nDims="2"
NVT="0"
qsub -v den=${den},eps=${eps},lambda=${lambda},restart=${restart},Lx=${Lx},Ly=${Ly},Lz=${Lz},Cx=${Cx},Cy=${Cy},Cz=${Cz},nDims=${nDims},NVT=${NVT} submit_gc_pillar.sh
