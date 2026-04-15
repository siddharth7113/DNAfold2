#!/bin/bash
cd scoring
mv top1.pdb CG.pdb
cp CG.pdb ../rebuild
cd ..
cd rebuild
bash run.sh
mv All_atom.pdb All_atom_top1.pdb
mv CG.pdb CG_top1.pdb
mv sec_struc.dat sec_struc_top1.dat
mv CG_top1.pdb ../../result/CG_strucure
mv All_atom_top1.pdb ../../result/All_atom_structure
mv sec_struc_top1.dat ../../result/Secondary_structure
cd ..
