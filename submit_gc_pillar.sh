#!/bin/sh
#$ -N GCpil_16.4
# Send mail to me:
#$ -M tretyakov@mpip-mainz.mpg.de
#$ -q Single.q
# Mail at beginning/end/on suspension:
#$ -m be
# The job is located in the current working directory:
#$ -cwd

#rm -rf /data/isilon/tretyakov/mf_lb/L${Lx}-${Ly}-${Lz}_C${Cx}/
#mkdir /data/isilon/tretyakov/mf_lb/L${Lx}-${Ly}-${Lz}_C${Cx}/
./gen_2015_03_13.out -path pckr160 -folder L${Lx}-${Ly}-${Lz}_C${Cx}/ \
-den ${den} \
-restart 1 \
-lambda ${lambda} \
-Lbox ${Lx} ${Ly} ${Lz} \
-corr ${Cx} ${Cy} ${Cz} \
-eps ${eps} \
-nMol 2000 \
-nDims ${nDims} \
-NVT ${NVT}
