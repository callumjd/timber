#!/bin/bash

##module purge
##module load Conda
##conda activate /home/dicksca3/.conda/envs/sstmap 

run_hsa -i ../prot_UNL.prmtop -t ../run_2/wat_image.nc -l ../UNL.pdb -f 40000 -o sstmap 

cd ./SSTMap_HSA

~/scripts/wat_hsa_output.py -i sstmap_hsa_summary.txt

