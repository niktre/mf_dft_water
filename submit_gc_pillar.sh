#!/bin/sh
#$ -N GCpil_13.4
# Send mail to me:
#$ -M tretyakov@mpip-mainz.mpg.de
#$ -q Single.q
# Mail at beginning/end/on suspension:
#$ -m be
# The job is located in the current working directory:
#$ -cwd

#rm -rf /data/isilon/tretyakov/mf_lb/e${eps}_den${den}_L${Lx}-${Ly}-${Lz}_C${Cx}_D${Dx}/
#mkdir /data/isilon/tretyakov/mf_lb/e${eps}_den${den}_L${Lx}-${Ly}-${Lz}_C${Cx}_D${Dx}/
./gen_2014_11_26.out -path pckr160 -folder e${eps}_den${den}_L${Lx}-${Ly}-${Lz}_C${Cx}_D${Dx}/ \
-den ${den} \
-restart 1 \
-lambda ${lambda} \
-Lbox ${Lx} ${Ly} ${Lz} \
-corr ${Cx} ${Cy} ${Cz} \
-cav ${Dx} ${Dy} ${Dz} \
-eps ${eps} \
-nMol 2000 \
-nDims ${nDims} \
-NVT ${NVT}
