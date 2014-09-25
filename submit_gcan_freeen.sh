#!/bin/sh
#$ -N GCpil_10.6
# Send mail to me:
#$ -M tretyakov@mpip-mainz.mpg.de
#$ -q Single.q
# Mail at beginning/end/on suspension:
#$ -m be
# The job is located in the current working directory:
#$ -cwd

mkdir /data/isilon/tretyakov/mf_lb/e${eps}_den${den}_L${Lx}-${Ly}-${Lz}_N${Nx}-${Ny}-${Nz}_C${Cx}_D${Dx}/
./gcan_2014_07_08.out -path pckr160 -folder e${eps}_den${den}_L${Lx}-${Ly}-${Lz}_N${Nx}-${Ny}-${Nz}_C${Cx}_D${Dx}/ \
-den ${den} \
-restart 0 \
-lambda 0.0010 \
-Lbox ${Lx} ${Ly} ${Lz} \
-grid ${Nx} ${Ny} ${Nz} \
-corr ${Cx} ${Cy} ${Cz} \
-cav ${Dx} ${Dy} ${Dz} \
-eps ${eps} \
-nMol 2000 \
-nDims ${nDims} \
-NVT ${NVT}
