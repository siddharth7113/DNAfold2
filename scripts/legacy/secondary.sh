#!/bin/bash
#SBATCH -n 16
rm -rf secondary
mkdir secondary
cp secondary.c sx.c cat.c all.sh op secondary
cd scoring
cp top1.pdb ch.dat ../secondary
cd ..
cp ch.c sdCG_ch.c op.c secondary
cd secondary
gcc -Wall secondary.c -o secondary -lm
./secondary
cp DNA_type ../
cp state.dat ../model
g++ ch.c
./a.out
gcc -O3 -Wall -fopenmp op.c -o op -lm
./op
bash all.sh
cd ..
