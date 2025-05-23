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

def xtb_torsion_scan(rd_mol,torsion_list,n_steps,step_size,gbsa_flag,output_name):

    # Ligand charge
    formal_charge=int(Chem.GetFormalCharge(rd_mol))

    # Write XYZ file
    with open('tmpAni.xyz','w') as f_out:
        f_out.write('%d\n' % (len(rd_mol.GetAtoms())))
        f_out.write('comment\n')
        for i in range(0,len(rd_mol.GetAtoms())):
            pos=rd_mol.GetConformer().GetAtomPosition(i)
            f_out.write('%s %lf %lf %lf\n' % (rd_mol.GetAtomWithIdx(i).GetSymbol(),pos.x,pos.y,pos.z))

    ene_values=[]
    dihed_values=[]
    sd_mols=[]

    start_angle=rdMolTransforms.GetDihedralDeg(rd_mol.GetConformer(),torsion_list[0]-1,torsion_list[1]-1,torsion_list[2]-1,torsion_list[3]-1)

    # xtb control file
    with open('cntrl_xtb.in','w') as f:
        f.write('$constrain\n')
        f.write('force constant=5\n')
        f.write('dihedral: %d,%d,%d,%d,auto\n' % (torsion_list[0],torsion_list[1],torsion_list[2],torsion_list[3]))
        f.write('$scan\n')
        f.write('1: %lf,%lf,%d\n' % (start_angle,start_angle+(n_steps*step_size),n_steps))

    if gbsa_flag==None:
        os.system('xtb tmpAni.xyz --opt -c %d -I cntrl_xtb.in > xtb_raw_out' % (formal_charge))
    else:
        os.system('xtb tmpAni.xyz --opt -c %d -I cntrl_xtb.in --alpb %s > xtb_raw_out' % (formal_charge,gbsa_flag))

    # get the energies and convert to kcal/mol
    with open('xtbscan.log','r') as f:
        for line in f:
            if 'energy:' in line:
                ene_values.append(float(line.split()[1]))
            elif 'SCF done' in line:
                ene_values.append(float(line.split()[2]))

    zero_point=min(ene_values)
    for i in range(0,len(ene_values)):
        ene_values[i]=(ene_values[i]-zero_point)*627.503

    # convert xyz structures to sdf
    with open('xtbscan.log','r') as f:
        data=f.readlines()

    for i in range(0,n_steps):
        local_mol=Chem.Mol(rd_mol)

        local_data=data[(i*len(rd_mol.GetAtoms()))+2+(i*2):(i*len(rd_mol.GetAtoms()))+2+(i*2)+len(rd_mol.GetAtoms())]

        for j in range(0,len(rd_mol.GetAtoms())):
            x=float(local_data[j].split()[1])
            y=float(local_data[j].split()[2])
            z=float(local_data[j].split()[3])

            local_mol.GetConformer().SetAtomPosition(j,(x,y,z))

        sd_mols.append(local_mol)
        dihed_values.append(rdMolTransforms.GetDihedralDeg(local_mol.GetConformer(),torsion_list[0]-1,torsion_list[1]-1,torsion_list[2]-1,torsion_list[3]-1))

        writer=SDWriter(output_name+'.sdf')
        counter=0
        for mol in sd_mols:
            this_scan=('%d-%d-%d-%d' % (torsion_list[0],torsion_list[1],torsion_list[2],torsion_list[3]))

            pm=PropertyMol(mol)
            pm.SetProp('xtb-kcal',ene_values[counter])
            pm.SetProp('torsion-angle',dihed_values[counter])
            pm.SetProp('torsion-scanned',this_scan)
            writer.write(pm)
            counter+=1
        writer.flush()

    # output the energy profile
    with open(output_name+'.dat','w') as f_out:
        f_out.write('# scan results %d %d %d %d\n' % (torsion_list[0],torsion_list[1],torsion_list[2],torsion_list[3]))
        for i in range(0,len(dihed_values)):
            f_out.write('%lf %lf\n' % (float(dihed_values[i]),float(ene_values[i])))

    os.system('rm tmpAni.xyz cntrl_xtb.in charges wbo xtbopt.log xtbscan.log xtbopt.xyz xtbrestart xtb_raw_out .xtboptok xtbtopo.mol')

###############################################################################
def run(input_sdf,output_name,torsion_index=None,torsion_smarts=None,gbsa=None):

    # Error checks
    if not os.path.isfile(input_sdf):
        print('Error: cannot find %s!\n' % (input_sdf))
        sys.exit()

    if torsion_index and torsion_smarts:
        print('Error: supply only torsion indices or smarts!\n')
        sys.exit()

    # Load SDF
    mol=Chem.SDMolSupplier(input_sdf,removeHs=False)[0]

    ######################
    # Submit torsion scan
    ######################
    if torsion_smarts and not torsion_index:
        torsion=get_torsion_smarts(mol,torsion_smarts)
        torsion=[x+1 for x in torsion]
    elif torsion_index and not torsion_smarts:
        torsion=torsion_index

    xtb_torsion_scan(mol,torsion,36,10,gbsa,output_name)

###############################################################################
## Parse arguments ##
###############################################################################

def main(argv=None):
    ## Command line arguments ##
    parser = argparse.ArgumentParser(description='XTB dihedral scan.\n')
    parser.add_argument('-i',help='Input SDF mol file',required=True)
    parser.add_argument('-n',help='Output name',required=False,default='torsion_xtb')
    parser.add_argument('-s',help='Scan dihedral flag',action='store',nargs=4,required=False)
    parser.add_argument('-smarts',help='Smarts pattern torsion match',action='store',required=False)
    parser.add_argument('-gbsa',help='Solvent model',action='store',required=False)
    args=vars(parser.parse_args())

###############################################################################
## Run selected protocols ##
###############################################################################
    run(input_sdf=args['i'],output_name=args['n'],torsion_index=args['s'],torsion_smarts=args['smarts'],gbsa=args['gbsa'])

# MAIN
if __name__ == '__main__':
    main()
    
