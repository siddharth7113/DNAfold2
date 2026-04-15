#!/bin/bash
#SBATCH -n 10
rm -rf model
mkdir model
bash initial.sh
cp TiRNA_remc config1.dat model
cd model
./TiRNA_remc
cd ..
cp t1.c model
cd model
g++ t1.c
./a.out
cd ..
bash scoring.sh
bash secondary.sh
bash optimize.sh
bash re_scoring.sh
bash rebuild.sh
g++ tm.c
./a.out
