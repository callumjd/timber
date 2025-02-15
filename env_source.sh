#!/bin/bash

module purge

module load uge
module load Amber/22-AmberTools22-CUDA11
module load OpenEye
module load PythonDS
conda activate /usr/prog/cadd/amber_tools/alchemistry2

export PYTHONPATH=$PYTHONPATH:/usr/prog/cadd/amber_tools/timber/versions/0.1/

export PATH=$PATH:/usr/prog/cadd/amber_tools/timber/versions/0.1/bin/
