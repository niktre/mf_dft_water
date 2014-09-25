#!/bin/bash
#Lx="16.4"
Lx="10.6"
Ly=${Lx}
Lz="16.4"
#Nx="164"
Nx="106"
Ny=${Nx}
Nz="164"
Cx="7.0"
Dx="3.6"
Cy=${Cx}
Dy=${Dx}
Cz="2."
Dz="0."
den="33.40"
nDims="2"
NVT="0"

#for eps in "0.01" "0.02" "0.03" "0.04" "0.05" "0.06" "0.07" "0.08" "0.09" "0.10" "0.11" "0.12" "0.13" "0.15" "0.18" "0.21" "0.24"
#for eps in "0.09" "0.10" "0.11" "0.12" "0.13" "0.15"
for eps in "0.15"
#for eps in "0.02" "0.03" "0.04" "0.05" "0.06" "0.07" "0.08" "0.09" "0.10" "0.11" "0.12" "0.13" "0.15" "0.18" "0.21" "0.24"
do
	qsub -v den=${den},eps=${eps},Lx=${Lx},Ly=${Ly},Lz=${Lz},Nx=${Nx},Ny=${Ny},Nz=${Nz},Cx=${Cx},Cy=${Cy},Cz=${Cz},Dx=${Dx},Dy=${Dy},Dz=${Dz},nDims=${nDims},NVT=${NVT} submit_gcan_freeen.sh
done
