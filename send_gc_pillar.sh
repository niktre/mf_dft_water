#!/bin/bash
Lx=$1
Ly=${Lx}
Lz="16.4"
Cx=$2
Cy=${Cx}
Cz="2."
den="33.40"
lambda=$3
nDims="2"
NVT="0"
eps="0.15"
qsub -v den=${den},eps=${eps},lambda=${lambda},Lx=${Lx},Ly=${Ly},Lz=${Lz},Cx=${Cx},Cy=${Cy},Cz=${Cz},nDims=${nDims},NVT=${NVT} submit_gc_pillar.sh
