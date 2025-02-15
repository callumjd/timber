#!/usr/bin/env python
import argparse
import glob
import sys
import os
import pandas as pd
import numpy as np

#
# Callum Dickson, NIBR CADD: callum.dickson@novartis.com
#

###############################################################################
## Routines ##
###############################################################################

def check_file(file_in):
    output=False
    if os.path.exists(file_in) and os.path.getsize(file_in)>0:
        output=True
    return output

def mapping_tuples(csv_file):

    map_list=[]

    # assumes first line are comments
    with open(csv_file,'r') as f:
        f.readline()
        for line in f:
            if line[0]!='#':
                lig1=line.split(',')[0].strip()
                if len(line.split(','))==2:
                    lig2=line.split(',')[1].strip()
                    map_list.append((lig1,lig2))
                else:
                    map_list.append((lig1))

    return map_list

def write_cluster(file_name,mode,lig,ff,sample,ensemble,gbsa):

    text1='''
#!/bin/bash
#$ -N auto-parmfit
#$ -j y
#$ -cwd
#$ -S /bin/bash
#$ -l h_rt=259100,m_mem_free=2G 

module load PythonDS
module load Amber
conda activate /usr/prog/cadd/amber_tools/alchemistry

'''

    with open(file_name,'w') as f:
        f.write(text1)
        f.write('sed -i \"s/UNL/%s/g\" make_lig.leap\n' % (lig))
        f.write('/usr/prog/cadd/amber_tools/alchemistry/bin/tleap -f make_lig.leap>out\n')
        if mode=='fragment':
            if gbsa:
                f.write('python /usr/prog/cadd/amber_tools/scripts/fragment_torfit.py -i %s.sdf -ff %s -gbsa\n' % (lig,ff))
            else:
                f.write('python /usr/prog/cadd/amber_tools/scripts/fragment_torfit.py -i %s.sdf -ff %s\n' % (lig,ff))

        elif mode=='torfit':
            if gbsa:
                if sample:
                    f.write('python /usr/prog/cadd/amber_tools/scripts/torfit.py -i %s.sdf -p prmtop -sample -gbsa -abc\n' % (lig))
                elif ensemble:
                    f.write('python /usr/prog/cadd/amber_tools/scripts/torfit.py -i %s.sdf -p prmtop -ensemble -gbsa -abc\n' % (lig))
                else:
                    f.write('python /usr/prog/cadd/amber_tools/scripts/torfit.py -i %s.sdf -p prmtop -gbsa -abc\n' % (lig))
            else:
                if sample:
                    f.write('python /usr/prog/cadd/amber_tools/scripts/torfit.py -i %s.sdf -p prmtop -sample -abc\n' % (lig))
                elif ensemble:
                    f.write('python /usr/prog/cadd/amber_tools/scripts/torfit.py -i %s.sdf -p prmtop -ensemble -abc\n' % (lig))
                else:
                    f.write('python /usr/prog/cadd/amber_tools/scripts/torfit.py -i %s.sdf -p prmtop -abc\n' % (lig))

        f.write('\n')

###############################################################################
# RUN
###############################################################################

def run(input_csv,mode,ff,sample,ensemble,gbsa):

    if not check_file(input_csv):
        print('Error: cannot find %s!\n' % (input_csv))
        sys.exit()

    map_list=mapping_tuples(input_csv)

    lig_name={'core':'LIG','sec_lig':'MOD'}

    for pair in map_list:
        if isinstance(pair,tuple):
            dir_name=pair[0]+'~'+pair[1]
            lig_list=['core','sec_lig']
        elif isinstance(pair,str): 
            dir_name=pair
            lig_list=['core']
        else:
            print('Cannot process map file!\n')
            sys.exit()

        # submit the jobs
        os.chdir(dir_name)

        for lig in lig_list:
            os.chdir(lig)

            write_cluster('parmfit.x',mode,lig_name[lig],ff,sample,ensemble,gbsa)
            os.system('qsub parmfit.x')

            os.chdir('../')

        os.chdir('../')

###############################################################################
## Parse arguments ##
###############################################################################

def main(argv=None):
    ## Command line arguments ##
    parser = argparse.ArgumentParser(description='Submit automated small molecule parameter fitting\n')
    
    parser.add_argument('-i','--input',dest='i',help='CSV file with ligand mappings',type=str,required=True)
    parser.add_argument('-m','--mode',dest='m',help='Mode of execution (default: fragment)',choices=['fragment','torfit'],default='fragment',required=False)
    parser.add_argument('-ff',dest='ff',help='Small molecule force field (default: gaff2)',choices=['gaff','gaff2','frosst'],default='gaff2',required=False)
    parser.add_argument('-sample', help='Flag to sample one extended conformer (torfit)',required=False,default=False,action='store_true')
    parser.add_argument('-ensemble', help='Flag to sample multiple extended conformers (torfit)',required=False,default=False,action='store_true')
    parser.add_argument('-gbsa',dest='gbsa',help='Use implicit water model',action='store_true',required=False,default=False)

    args=vars(parser.parse_args())

###############################################################################
## Run selected protocols ##
###############################################################################

    run(input_csv=args['i'],mode=args['m'],ff=args['ff'],sample=args['sample'],ensemble=args['ensemble'],gbsa=args['gbsa'])

# MAIN
if __name__ == '__main__':
    main()

