#!/usr/bin/env python
import sys
import os
import numpy as np
import math
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops
from rdkit.Chem import rdMolTransforms
import argparse

#
# Callum Dickson, NIBR CADD: callum.dickson@novartis.com
#

# WARNING: Will only work for covalent CYS and ASP residues
# Other residues are possible but will need to be added to the script

##############################################################################
## Routines ##
###############################################################################

def check_file(file_in):
    output=False
    if os.path.exists(file_in) and os.path.getsize(file_in)>0:
        output=True
    return output

def check_complete(input_file):

    with open(input_file,'r') as f:
        data=f.readlines()
    if 'Normal termination of Gaussian' in data[-1]:
        return True
    else:
        return False

def make_g09(input_sdf):

    mol=Chem.SDMolSupplier(input_sdf,removeHs=False)[0]
    if not mol:
        print('Cannot load %s!\n' % (input_sdf))
        sys.exit()

    chg=int(rdmolops.GetFormalCharge(mol))

    with open('opt.com','w') as f:
        f.write('%nprocshared=8\n')
        f.write('%mem=32GB\n')
        f.write('%chk=opt.chk\n')
        # fixme: opt
        f.write('#P b3lyp/6-31G* opt\n')
        f.write('\n')
        f.write('opt\n')
        f.write('\n')
        f.write('%d 1\n' % (chg))

        for i in range(0,len(mol.GetAtoms())):
            pos=mol.GetConformer().GetAtomPosition(i)
            f.write(' %s    %lf %lf %lf\n' % (mol.GetAtomWithIdx(i).GetSymbol(),pos.x,pos.y,pos.z))
        f.write('\n')

    with open('charge.com','w') as f:
        f.write('%nprocshared=8\n')
        f.write('%mem=32GB\n')
        f.write('%chk=charge.chk\n')
        f.write('#P HF/6-31G* SCF=Tight Pop=MK IOp(6/33=2) geom=checkpoint\n')
        f.write('\n')
        f.write('charge\n')
        f.write('\n')
        f.write('%d 1\n' % (chg))
        f.write('\n')

    with open('qm.x','w') as f:
        f.write('#!/bin/bash\n')
        f.write('\n')
        f.write('#$ -N QM\n')
        f.write('#$ -j y\n')
        f.write('#$ -cwd\n')
        f.write('#$ -S /bin/bash\n')
        f.write('#$ -V\n')
        f.write('#$ -l h_rt=259100,m_mem_free=4G\n')
        f.write('#$ -pe smp 8\n')
        f.write('\n')
        f.write('module load Gaussian/09-foss-2018a\n')
        f.write('\n')
        f.write('g09 <opt.com>opt.log\n')
        f.write('\n')
        f.write('cp opt.chk charge.chk\n')
        f.write('\n')
        f.write('g09 <charge.com>charge.log\n')
        f.write('\n')

    os.system('qsub qm.x')

def resp_fit(input_sdf,charge_log):

    mol=Chem.SDMolSupplier(input_sdf,removeHs=False)[0]
    if not mol:
        print('Cannot load %s!\n' % (input_sdf))
        sys.exit()

    chg=int(rdmolops.GetFormalCharge(mol))
    n_atoms=int(len(mol.GetAtoms()))

    # get the esp.sh and readit.f files
    os.system('cp /usr/prog/cadd/amber_tools/extra-scripts/covAMBER/esp.sh .')
    os.system('cp /usr/prog/cadd/amber_tools/extra-scripts/covAMBER/replace.f .')

    nesp=None
    ctr=0
    with open(charge_log,'r') as f:
        for line in f:
            if 'Fit  ' in line:
                ctr+=1
    nesp=ctr
    if not nesp:
        print('Error: cannot count ESP fits!\n')
        sys.exit()

    os.system('cp replace.f readit.f')
    os.system('sed -i "s,NATOM,%d,g" readit.f' % (n_atoms))
    os.system('sed -i "s,NESP,%d,g" readit.f' % (nesp))

    # Write the resp input files
    CH2=mol.GetSubstructMatches(AllChem.MolFromSmarts('[C;H2]([H])[H]'))
    CH3=mol.GetSubstructMatches(AllChem.MolFromSmarts('[C;H3]([H])([H])[H]'))

    nme='[H][C;H3]([H])([H])[N;H1]([H])[C](=O)-[C;H1]([H])-[C;H2]([H])([H])-[S,#6]'
    ace='[H][C;H3]([H])([H])[C](=O)[N;H1]([H])-[C;H1]([H])-[C;H2]([H])([H])-[S,#6]'

    freeze_nme=list(mol.GetSubstructMatches(AllChem.MolFromSmarts(nme))[0][0:6])
    freeze_ace=list(mol.GetSubstructMatches(AllChem.MolFromSmarts(ace))[0][0:6])

    freeze_idx=freeze_nme+freeze_ace

    carbons=[]
    for hit in CH2:
        carbons.append(hit[0])
    for hit in CH3:
        carbons.append(hit[0])

    first_H=[]
    for hit in CH2:
        sub=sorted(hit[1:])
        carbons.append(sub[0])
    for hit in CH3:
        sub=sorted(hit[1:])
        carbons.append(sub[0])

    partner={}
    for hit in CH2:
        sub=sorted(hit[1:])
        partner.update({sub[1]:(sub[0])})
    for hit in CH3:
        sub=sorted(hit[1:])
        partner.update({sub[1]:(sub[0])})
        partner.update({sub[2]:(sub[0])})

    # File 1.
    # Any symmetry will not be reflected (correct ..) 
    with open('fit_resp.in','w') as f:
        f.write('resp run #1\n')
        f.write(' &cntrl\n')
        f.write(' ihfree=1,\n')
        f.write(' qwt=0.0005,\n')
        f.write(' iqopt=2,\n')
        f.write(' /\n')
        f.write('\n')
        f.write('mol\n')
        f.write('    %d%5d\n' % (chg,n_atoms))
        for i in range(0,n_atoms):
            atm_num=int(mol.GetAtomWithIdx(i).GetAtomicNum())
            if i in freeze_idx:
                f.write('%5d   -1\n' % (atm_num))
            else:
                f.write('%5d    0\n' % (atm_num))
        f.write('\n')

    # File 2.
    with open('fit_resp2.in','w') as f:
        f.write('resp run #2\n')
        f.write(' &cntrl\n')
        f.write(' ihfree=1,\n')
        f.write(' qwt=0.001,\n')
        f.write(' iqopt=2,\n')
        f.write(' /\n')
        f.write('\n')
        f.write('mol\n')
        f.write('    %d%5d\n' % (chg,n_atoms))
        for i in range(0,len(mol.GetAtoms())):
            atm_num=int(mol.GetAtomWithIdx(i).GetAtomicNum())
            if i in freeze_idx:
                f.write('%5d   -1\n' % (atm_num))
            elif i in carbons:
                f.write('%5d    0\n' % (atm_num))
            elif i in first_H:
                f.write('%5d    0\n' % (atm_num))
            elif i in partner:
                f.write('%5d   %d\n' % (atm_num,partner[i]+1))
            else:
                f.write('%5d   -1\n' % (atm_num))
        f.write('\n')

    nme_qin=[0.097600,-0.149000,0.097600,0.097600,-0.415700,0.271900]
    ace_qin=[0.112300,-0.366200,0.112300,0.112300,0.597200,-0.567900]

    nme_dict={}
    ctr=0
    for val in freeze_nme:
        nme_dict.update({val:nme_qin[ctr]})
        ctr+=1

    ace_dict={}
    ctr=0
    for val in freeze_ace:
        ace_dict.update({val:ace_qin[ctr]})
        ctr+=1

    all_charges=[]
    for i in range(0,n_atoms):
        if i in freeze_nme:
            all_charges.append(nme_dict[i])
        elif i in freeze_ace:
            all_charges.append(ace_dict[i])
        else:
            all_charges.append(0.000000)

    # Write the qin file
    with open('qin.qin','w') as f:
    
        for i in range(0,int(math.floor(n_atoms/8))):

            line=''
            for val in all_charges[(i*8):(i*8)+8]:
                line=line+('%10.6f' % (val))
            line=line+'\n'

            f.write(line)

        line=''
        for val in all_charges[(i*8)+8:]:
            line=line+('%10.6f' % (val))
        line=line+'\n'

        f.write(line)

    if check_file('esp.dat'):
        os.system('rm esp.dat')

    os.system('bash esp.sh charge.log')
    os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/resp -O -i fit_resp.in -o fit_resp.out -p fit_resp.pch -t fit_resp.chg -e esp.dat -q qin.qin')
    os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/resp -O -i fit_resp2.in -o fit_resp2.out -p fit_resp2.pch -t fit_resp2.chg -q fit_resp.chg -e esp.dat')

def build_final(input_sdf,start_pdb,start_off,start_chg):

    mol=Chem.SDMolSupplier(input_sdf,removeHs=False)[0]
    if not mol:
        print('Cannot load %s!\n' % (input_sdf))
        sys.exit()

    chg=int(rdmolops.GetFormalCharge(mol))
    n_atoms=int(len(mol.GetAtoms()))

    nme='[H][C;H3]([H])([H])[N;H1]([H])[C](=O)-[C;H1]([H])-[C;H2]([H])([H])-[S,#6]'
    ace='[H][C;H3]([H])([H])[C](=O)[N;H1]([H])-[C;H1]([H])-[C;H2]([H])([H])-[S,#6]'
    # Look for cysteine or aspartic acid
    sulf='[H][N;H1]C([H])([#6]=O)C([H])([H])[#16]'
    asp='[H][N;H1]C([H])([#6]=O)C([H])([H])[#6](=O)[O;X2]'

    drop_nme=list(mol.GetSubstructMatches(AllChem.MolFromSmarts(nme))[0][0:6])
    drop_ace=list(mol.GetSubstructMatches(AllChem.MolFromSmarts(ace))[0][0:6])
    drop_all=drop_nme+drop_ace

    just_co=list(mol.GetSubstructMatches(AllChem.MolFromSmarts(nme))[0])[6]
    just_nh=list(mol.GetSubstructMatches(AllChem.MolFromSmarts(ace))[0])[6]

    co_idx=None
    nh_idx=None
    ctr=0
    for i in range(0,n_atoms):
        if i not in drop_all:
            if i==just_co:
                co_idx=ctr
            elif i==just_nh:
                nh_idx=ctr
            ctr+=1

    if mol.HasSubstructMatch(AllChem.MolFromSmarts(sulf)):
        keep_sulf=list(mol.GetSubstructMatches(AllChem.MolFromSmarts(sulf))[0])
        sulf_types=['H','N','CX','H1','C','O','2C','H1','H1','S']
    elif mol.HasSubstructMatch(AllChem.MolFromSmarts(asp)):
        keep_sulf=list(mol.GetSubstructMatches(AllChem.MolFromSmarts(asp))[0])
        sulf_types=['H','N','CX','H1','C','O','2C','H1','H1','CO','O2','OS']
    else:
        print('Error: cannot parse covalent attachment!\n')
        sys.exit()

    sulf_dict={}
    ctr=0
    for val in keep_sulf:
        sulf_dict.update({val:sulf_types[ctr]})
        ctr+=1

    # write a COV.pdb file, skipping drop_nme and drop_ace atoms
    with open(start_pdb,'r') as f:
        pdb=f.readlines()

    atom_names=[]
    with open(start_pdb,'r') as f:
        for line in f:
            name=line.split()[2]
            atom_names.append(name)

    with open('COV.pdb','w') as f:
        for i in range(0,len(pdb)):
            if i not in drop_all:
                f.write(pdb[i].replace('UNL','COV'))

    # save all types; update CYS region
    at_types=[]
    with open(start_off,'r') as f:
        f.readline()
        f.readline()
        f.readline()
        ctr=0
        for line in f:
            if ctr<n_atoms:
                at=line.split()[1].strip('"').strip()
                at_types.append(at)
                ctr+=1

    for k,v in sulf_dict.items():
        at_types[k]=v

    # get all charges
    at_charges=[]
    with open(start_chg,'r') as f:
        for line in f:
            for val in line.split():
                at_charges.append(float(val.strip()))

    # write a file for leap
    with open('make_final.leap','w') as f:
        f.write('COV=loadpdb COV.pdb\n')
        
        for i in range(0,n_atoms):
            if i not in drop_all:
                f.write('set COV.1.%s type "%s"\n' % (atom_names[i],at_types[i]))

        for i in range(0,n_atoms):
            if i not in drop_all:
                f.write('set COV.1.%s charge %lf\n' % (atom_names[i],at_charges[i]))

        for bond in mol.GetBonds():
            a=int(bond.GetBeginAtomIdx())
            b=int(bond.GetEndAtomIdx())

            if (a not in drop_all) and (b not in drop_all):
                f.write('bond COV.1.%s COV.1.%s\n' % (atom_names[a],atom_names[b]))

        f.write('saveoff COV COV.off\n')
        f.write('quit\n')

    os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/tleap -f make_final.leap>leap.out')

    with open('leap.out','r') as f:
        data=f.readlines()

    if 'Exiting LEaP: Errors = 0' in data[-1]:
        print('Build was successful\n')
    else:
        print('Error: leap build did not complete!\n')

    # fill the connect array
    if 'Exiting LEaP: Errors = 0' in data[-1]:
        with open('COV.off','r') as f:
            off=f.readlines()

        ctr=0
        split_line=None
        for line in off:
            if 'unit.connect array int' in line:
                split_line=ctr
            ctr+=1

        with open('COV.off','w') as f:
            for i in range(0,split_line+1):
                f.write(off[i])
            f.write(' %d\n' % (nh_idx+1))
            f.write(' %d\n' % (co_idx+1))
            for i in range(split_line+3,len(off)):
                f.write(off[i])

def run(input_sdf):

    if not check_file(input_sdf):
        print('\nError: cannot find %s!\n' % (input_sdf))
        sys.exit()

    # Prepare and submit the QM
    if not check_file('opt.com'):
        os.system('/usr/prog/cadd/amber_tools/scripts/run_ligprep.py -i %s -ff frosst' % (input_sdf))

        if not check_file('UNL.sdf'):
            print('\nError: initial parameter generation failed!\n')
            sys.exit()

        print('\nPrepare and submit QM\n')
        make_g09('UNL.sdf')

    elif check_file('opt.com') and not check_file('charge.log'):
        print('\nError: is QM still running?\n')
        sys.exit()

    # QM is complete: process
    if check_file('opt.com') and check_file('charge.log') and not check_file('fit_resp2.out'):
        if not check_complete('charge.log'):
            print('\nError: QM died!\n')
            sys.exit() 

        if check_complete('charge.log'):
            print('\nProcess charge.log\n')
            resp_fit('UNL.sdf','charge.log')

    # resp finished: build the final parameters
    if check_file('opt.com') and check_file('charge.log') and check_file('fit_resp2.out'):
        print('\nBuild final COV parameters: COV.pdb and COV.off\n')
        build_final('UNL.sdf','UNL.pdb','UNL.off','fit_resp2.chg')


###############################################################################
## Parse arguments ##
###############################################################################

def main(argv=None):
    ## Command line arguments ##
    parser = argparse.ArgumentParser(description='Submit short MD refinement with AMBER\n')
    parser.add_argument('-i','--input',dest='i',help='Input SDF of capped ligand (ACE,NME)',type=str,required=True)

    args=vars(parser.parse_args())

###############################################################################
## Run selected protocols ##
###############################################################################
    run(input_sdf=args['i'])

# MAIN
if __name__ == '__main__':
    main()

