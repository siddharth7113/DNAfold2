#!/bin/bash
cd ..
cd model
cp conf_0.dat Energy_0.dat ch.dat ../scoring
cd ..
cd scoring
rm -rf Cluster*
python Umin.py
gcc -Wall A_state.c -o A_state -lm
./A_state
g++ tc.c
./a.out
python cluster1.py
python cluster2.py
rm temp*.pdb
cp cluster3.py Cluster_0
cd Cluster_0
python cluster3.py
cp top1.pdb ..
