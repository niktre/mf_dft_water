#!/bin/bash
Lx=$1
Ly=${Lx}
Lz="16.4"
Cx=$2
Dx=$3
Cy=${Cx}
Dy=${Dx}
Cz="2."
Dz="0."
den="33.40"
lambda=$4
nDims="2"
NVT="0"
eps="0.15"
#Dx=$(( ${Lx}-${Cx} ))
qsub -v den=${den},eps=${eps},lambda=${lambda},Lx=${Lx},Ly=${Ly},Lz=${Lz},Cx=${Cx},Cy=${Cy},Cz=${Cz},Dx=${Dx},Dy=${Dy},Dz=${Dz},nDims=${nDims},NVT=${NVT} submit_gc_pillar.sh
