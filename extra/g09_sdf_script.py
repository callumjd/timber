#!/usr/bin/env python
import os
import sys
import numpy as np
import linecache
import argparse
from rdkit import Chem
from rdkit.Chem import SDWriter
from rdkit.Chem import AllChem
from rdkit.Chem.PropertyMol import *
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import rdmolops

###############################################################################
## Routines ##
###############################################################################

def check_complete(input_file):
    with open(input_file,'r') as f:
        data=f.readlines()
    if 'Normal termination of Gaussian 09' in data[-1]:
        return True
    else:
        return False

def extract_ene(input_file):
    i=0
    complete=[]
    with open(input_file,'r') as f:
        for line in f:
            if 'Optimization completed' in line:
                complete.append(i)
            i+=1

    ene_final=[]
    for val in range(0,len(complete)):
        j=0
        while True:
            line=linecache.getline(input_file,complete[val]-j)
            j+=1
            if 'SCF Done' in line:
                ene_final.append(float(line.split()[4]))
                break

    return ene_final

def zero_convert(dft_ene):
    kcal_out=[]
    min_val=min(dft_ene)

    for val in dft_ene:
        kcal_out.append((val+abs(min_val))*627.5095602)

    return kcal_out

def find_tor_idx(input_file):
    dihe_line=0
    with open(input_file,'r') as f:
        for line in f:
            dihe_line+=1
            if 'The following ModRedundant input section' in line:
                break

    dihe=linecache.getline(input_file,dihe_line+1)
    a=int(dihe.split()[1])
    b=int(dihe.split()[2])
    c=int(dihe.split()[3])
    d=int(dihe.split()[4])

    return [a,b,c,d]

def compile_struc(input_file,rd_sdf):

    i=0
    complete=[]
    with open(input_file,'r') as f:
        for line in f:
            if 'Optimization completed' in line:
                complete.append(i)
            i+=1

    struc=[]
    for val in range(0,len(complete)):
        j=0
        while True:
            line=linecache.getline(input_file,complete[val]-j)
            j+=1
            if 'Standard orientation' in line:
                struc.append(complete[val]-j)
                break

    out_mols=[]
    for val in range(0,len(struc)):
        local_mol=PropertyMol(rd_sdf)
        for j in range(0,len(rd_sdf.GetAtoms())):
            line=linecache.getline(input_file,struc[val]+6+j)
            x=float(line.split()[3])
            y=float(line.split()[4])
            z=float(line.split()[5])

            local_mol.GetConformer().SetAtomPosition(j,(x,y,z))

        out_mols.append(local_mol)

    return out_mols

def append_info(mol_list,tor_idx,ene_vals):

    this_scan=('%d-%d-%d-%d' % (tor_idx[0],tor_idx[1],tor_idx[2],tor_idx[3]))

    final_mols=[]
    counter=0
    for mol in mol_list:
        torsion_angle=rdMolTransforms.GetDihedralDeg(mol.GetConformer(),tor_idx[0]-1,tor_idx[1]-1,tor_idx[2]-1,tor_idx[3]-1)
        mol.SetProp('torsion-scanned',this_scan)
        mol.SetProp('torsion-angle',torsion_angle)
        mol.SetProp('kcal_energy',ene_vals[counter])

        final_mols.append(mol)

        counter+=1

    return final_mols

def write_jobscript(name,g09_name):

    text='''
#!/bin/bash

#$ -N QM 
#$ -j y
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -l h_rt=259100,m_mem_free=4G # note that total memory is 2 * 2G
#$ -binding linear:4
#$ -pe smp 4

#module load module load Gaussian/09-foss-2018a

'''

    with open(name,'w') as f:
        f.write(text)
        f.write('g09 <%s.com>%s.log\n' % (g09_name,g09_name))
        f.write('\n')

def write_g09(mol,name,torsion,gbsa=None):

    # Set keyworkds, mem and nproc
    # Hard coded here inside function
    mem='16GB'
    nprocshared=4
    level='opt=modredundant b3lyp/6-31g* EmpiricalDispersion=GD3'

    if gbsa:
        level=level+' scrf=(iefpcm,solvent=%s)' % (gbsa)

    # Torsion
    a=int(torsion[0])
    b=int(torsion[1])
    c=int(torsion[2])
    d=int(torsion[3])

    # Ligand charge
    mol_charge=int(Chem.GetFormalCharge(mol))

    with open(name,'w') as f:
        f.write('%%nprocshared=%d\n' % (nprocshared))
        f.write('%%mem=%s\n' % (mem))
        f.write('%%chk=%s.chk\n' % (name.split('.')[0]))
        f.write('# %s \n' % (level))
        f.write('\n')
        f.write('%s\n' % (name.split('.')[0]))
        f.write('\n')
        f.write('%d 1\n' % (mol_charge))

        for atom in range(0,mol.GetNumAtoms()):
            symbol=mol.GetAtomWithIdx(atom).GetSymbol()
            pos=mol.GetConformer().GetAtomPosition(atom)

            f.write('%s %lf %lf %lf\n' % (symbol,pos.x,pos.y,pos.z))

        f.write('\n')
        f.write('D %d %d %d %d S 36 10.0\n' % (a,b,c,d))
        f.write('\n')

def get_torsion_smarts(mol,smarts):

    smarts_query=AllChem.MolFromSmarts(smarts)
    if not mol.HasSubstructMatch(smarts_query):
        print('Error: cannot find smarts match!\n')
        sys.exit()

    matches=mol.GetSubstructMatches(smarts_query)
    if len(matches)<5:
        return matches[0]
    else:
        print('Error: found >4 torsion matches!\n')
        sys.exit()

def write_sdf_scan(info_mols,output_name,tor_idx):

    # Write SDF File
    print('Output structures: %s.sdf\n' % (output_name))
    writer=SDWriter(output_name+'.sdf')
    for mol in info_mols:
        writer.write(mol)
    writer.flush()

    # Write .dat File
    print('Output profile: %s.dat\n' % (output_name))
    with open(output_name+'.dat','w') as f:
        f.write('# Scan results %d %d %d %d\n' % (tor_idx[0],tor_idx[1],tor_idx[2],tor_idx[3]))
        for mol in info_mols:
            f.write('%lf %lf\n' % (float(mol.GetProp('torsion-angle')),float(mol.GetProp('kcal_energy'))))

###############################################################################
def run(input_sdf,output_name,logfile=None,torsion_index=None,torsion_smarts=None,gbsa=None):

    # Error checks
    if not os.path.isfile(input_sdf):
        print('Error: cannot find %s!\n' % (input_sdf))
        sys.exit()

    if logfile:
        if not os.path.isfile(logfile):
            print('Error: cannot find %s!\n' % (logfile))
            sys.exit()

    if torsion_index and torsion_smarts:
        print('Error: supply only torsion indices or smarts!\n')
        sys.exit()

    # Load SDF
    mol=Chem.SDMolSupplier(input_sdf,removeHs=False)[0]

    ######################
    # Submit torsion scan
    ######################
    if any([torsion_index,torsion_smarts]):
        if torsion_smarts and not torsion_index:
            torsion=get_torsion_smarts(mol,torsion_smarts)
            torsion=[x+1 for x in torsion]
        elif torsion_index and not torsion_smarts:
            torsion=torsion_index

        # Write .com file
        if gbsa=='h2o':
            gbsa='water'

        write_g09(mol,'opt.com',torsion,gbsa)

        # Submit
        print('\nSubmit torsion scan\n')
        write_jobscript('run.x','opt')
        os.system('qsub run.x')

    ######################
    # Parse torsion scan
    ######################
    if not any([torsion_index,torsion_smarts]):
        if not logfile:
            print('Error: no torsion indices or logfile supplied!\n')
            sys.exit()

        if not check_complete(logfile):
            print('Error: g09 %s not complete!\n' % (logfile))
            sys.exit()

        print('\nParsing %s scan\n' % (logfile))

        dft_raw=extract_ene(logfile)
        kcal_raw=zero_convert(dft_raw)
        tor_idx=find_tor_idx(logfile)
        scan_mols=compile_struc(logfile,mol)
        info_mols=append_info(scan_mols,tor_idx,kcal_raw)

        write_sdf_scan(info_mols,output_name,tor_idx)

###############################################################################
## Parse arguments ##
###############################################################################

def main(argv=None):
    ## Command line arguments ##
    parser = argparse.ArgumentParser(description='G09 dihedral scan.\n')
    parser.add_argument('-i',help='Input SDF mol file',required=True)
    parser.add_argument('-n',help='Output name',required=False,default='torsion_g09')
    parser.add_argument('-l',help='Gaussian log file to parse',required=False)
    parser.add_argument('-s',help='Scan dihedral flag',action='store',nargs=4,required=False)
    parser.add_argument('-smarts',help='Smarts pattern torsion match',action='store',required=False)
    parser.add_argument('-gbsa',help='Solvent model',action='store',required=False)
    args=vars(parser.parse_args())

###############################################################################
## Run selected protocols ##
###############################################################################
    run(input_sdf=args['i'],output_name=args['n'],logfile=args['l'],torsion_index=args['s'],torsion_smarts=args['smarts'],gbsa=args['gbsa'])

# MAIN
if __name__ == '__main__':
    main()
    
