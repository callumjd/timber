#!/usr/bin/env python

import argparse
import os
import sys
import glob as glob
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools
from rdkit.Chem import rdMolTransforms,rdmolops,SDWriter

###############################################################################

def setup_dir(all_sdf):

    for sdf in all_sdf:
        df=PandasTools.LoadSDF(sdf,removeHs=False)
        df=df.sort_values(by=['weight_water'],ascending=[False])
        df=df.reset_index(drop=True)

        name=str(df.loc[0,'ID'])

        if not os.path.exists(name):
            os.mkdir(name)
            os.chdir(name)

            PandasTools.WriteSDF(df,sdf,idName='ID',molColName='ROMol',properties=list(df.columns))

            for index,row in df.iterrows():
                dirname='conf_'+str(index)
                os.mkdir(dirname)
            
                writer=SDWriter(dirname+'/input.sdf')
                writer.write(row['ROMol'])
                writer.flush()

            os.chdir('../')

        else:
            print('Error: %s dir exists!\n' % (name))
            sys.exit()

def run_xtb(all_sdf):

    inp='''$write
   gbsa=true
$end

'''

    for sdf in all_sdf:
        suppl=Chem.SDMolSupplier(sdf,removeHs=False)
        name=suppl[0].GetProp('_Name')
        formal_chg=int(rdmolops.GetFormalCharge(suppl[0]))
 
        if os.path.isdir(name):
            os.chdir(name)
            
            for conf in glob.glob('conf_*'):
                os.chdir(conf)

                # XTB input
                with open('xtb.inp','w') as f:
                    f.write(inp)
    
                # XTB OPT
                os.system('xtb --gfn2 --opt vtight --acc 0.001 --chrg %d --alpb water --input xtb.inp --molden --namespace optimization input.sdf > optimization.output' % (formal_chg))
    
                # XTB QPLUS
                os.system('xtb --gfn2 --acc 0.001 --chrg %d --alpb water --namespace qplus optimization.xtbopt.sdf > qplus.output' % (formal_chg-1))

                # XTB FUKUI
                os.system('xtb --gfn2 --vfukui --acc 0.001 --chrg %d --alpb water --namespace vfukui optimization.xtbopt.sdf > vfukui.output' % (formal_chg))
    
                # VIPEA
                os.system('xtb --gfn2 --vipea --acc 0.001 --chrg %d --alpb water --namespace vipea optimization.xtbopt.sdf > vipea.output' % (formal_chg))
    
                # VOMEGA
                os.system('xtb --gfn2 --vomega --acc 0.001 --chrg %d --alpb water --namespace vomega optimization.xtbopt.sdf > vomega.output' % (formal_chg))

                os.chdir('../')

            os.chdir('../')

# A function to return the Boltzmann weighted average of values.
def Boltzmann_avg(values, weights):
    
    if len(values) != 0:            
        # Get a normalization factor for the new weights.
        norm = sum(weights)

        # Calculate the Boltzmann weighted average.
        B_avg = 0
        for x in range(len(values)):
            B_avg = B_avg + (values[x]*weights[x]/norm)

    # Return the final value.
    return B_avg

def gather_xtb_ene(weights):

    kcal_mol = 627.5

    all_conf=glob.glob('conf_*')
    all_conf=sorted(all_conf,key=lambda x: int(x.split('_')[1]))

    xtb_ene=[]
   
    try: 
        for c in all_conf:
            outfile_path = c+'/optimization.output'
            with open(outfile_path,'r') as f:
                for line in f:
                    if 'TOTAL ENERGY' in line:
                        E_total = float(line.split()[-3])*kcal_mol        
                        xtb_ene.append(E_total)

        if len(xtb_ene)==len(weights):
            return Boltzmann_avg(xtb_ene,weights)

    except:
        return 0

def gather_xtb_hlgap(weights):

    all_conf=glob.glob('conf_*')
    all_conf=sorted(all_conf,key=lambda x: int(x.split('_')[1]))

    hl=[]

    try:
        for c in all_conf:
            outfile_path = c+'/optimization.output'
            with open(outfile_path,'r') as f:
                for line in f:
                    if '          | HOMO-LUMO GAP' in line:
                        hl.append(float(line.split()[-3]))

        if len(hl)==len(weights):
            return Boltzmann_avg(hl,weights)

    except:
        return 0

def gather_xtb_dipole(weights):

    all_conf=glob.glob('conf_*')
    all_conf=sorted(all_conf,key=lambda x: int(x.split('_')[1]))

    dipole=[]

    try:
        for c in all_conf:
            outfile_path = c+'/optimization.output'
            outfile = open(outfile_path, 'r')
            outlines = outfile.readlines()
            outfile.close()  
    
            # Get the dipole moment.
            for i in range(len(outlines)):
                ol = outlines[i]
                if ol.startswith('molecular dipole:'):
                    dm = float(outlines[i+3].split()[-1])
                    dipole.append(dm)

        if len(dipole)==len(weights):
            return Boltzmann_avg(dipole,weights)

    except:
        return 0

# global electrophilicity index (eV)
def gather_xtb_ei(weights):

    all_conf=glob.glob('conf_*')
    all_conf=sorted(all_conf,key=lambda x: int(x.split('_')[1]))

    ei=[]

    try:
        for c in all_conf:
            outfile_path = c+'/vomega.output'
            with open(outfile_path,'r') as f:
                for line in f:
                    if 'Global electrophilicity index' in line:
                        ei.append(float(line.split()[-1]))

        if len(ei)==len(weights):
            return Boltzmann_avg(ei,weights)

    except:
        return 0

# electron affinity
def gather_xtb_ea(weights):

    all_conf=glob.glob('conf_*')
    all_conf=sorted(all_conf,key=lambda x: int(x.split('_')[1]))

    ea=[]

    try:
        for c in all_conf:
            outfile_path = c+'/vipea.output'
            with open(outfile_path,'r') as f:
                for line in f:
                    if 'delta SCC EA (eV)' in line:
                        ea.append(float(line.split()[-1]))

        if len(ea)==len(weights):
            return Boltzmann_avg(ea,weights)

    except:
        return 0

def gather_fukui_ensemble(atom_idx):

    all_conf=glob.glob('conf_*')
    all_conf=sorted(all_conf,key=lambda x: int(x.split('_')[1]))

    f=[]

    try:
        for c in all_conf:
            outfile_path = c+'/vfukui.output'
            fukui_file = open(outfile_path, 'r')
            fukui_lines = fukui_file.readlines()
            fukui_file.close()
    
            # Get the values.
            for j in range(len(fukui_lines)):
                if len(fukui_lines[j].split()) == 4 and fukui_lines[j].split()[1] == 'f(+)':
                    fukui = float(fukui_lines[j+1+atom_idx].split()[1])
                    f.append(fukui)

        return f

    except:
        return 0

def gather_charge_ensemble(atom_idx,basename):

    all_conf=glob.glob('conf_*')
    all_conf=sorted(all_conf,key=lambda x: int(x.split('_')[1]))

    q=[]

    try:
        for c in all_conf:
            outfile_path = c+'/'+basename+'.charges'
            charge_file = open(outfile_path, 'r')
            charge_lines = charge_file.readlines()
            charge_file.close()

            # Get the atomic charges.
            charges = []
            for c in charge_lines:
                charges.append(float(c.split()[0]))

            q.append(charges[atom_idx])

        return q

    except:
        return 0

def gather_sasa_ensemble(atom_idx):

    all_conf=glob.glob('conf_*')
    all_conf=sorted(all_conf,key=lambda x: int(x.split('_')[1]))

    sasa=[]

    try:
        for c in all_conf:
            outfile_path = c+'/optimization.output'
            sasa_file = open(outfile_path, 'r')
            sasa_lines = sasa_file.readlines()
            sasa_file.close()

            # Find the start of the SASA section.
            for j in range(len(sasa_lines)):
                if sasa_lines[j].startswith(' * generalized Born model'):
            
                    # Get the SASA value
                    s = float(sasa_lines[j+3+atom_idx].split()[4])
                    sasa.append(s)

        return sasa

    except:
        return 0

def gather_xtb(all_sdf,smarts):

    df=pd.DataFrame(columns=['Id','XTB_Energy','Homo-Lumo','Dipole','EI','EA'])

    # f+ Fukui indices
    for i in range(0,len(Chem.MolFromSmarts(smarts).GetAtoms())):
        df['f+_'+str(i)]=None

    # partial charges
    for i in range(0,len(Chem.MolFromSmarts(smarts).GetAtoms())):
        df['q_'+str(i)]=None

    # atomic charges with an extra electron
    for i in range(0,len(Chem.MolFromSmarts(smarts).GetAtoms())):
        df['q+_'+str(i)]=None

    # atomic SASA
    for i in range(0,len(Chem.MolFromSmarts(smarts).GetAtoms())):
        df['sasa_'+str(i)]=None

    counter=0
    for sdf in all_sdf:
        suppl=Chem.SDMolSupplier(sdf,removeHs=False)
        name=suppl[0].GetProp('_Name')

        if os.path.isdir(name):

            df.at[counter,'Id']=name

            os.chdir(name)

            # load the re-ordered SDF
            mp=PandasTools.LoadSDF(sdf)
            mp['weight_water']=mp['weight_water'].astype(float)
            weights=list(mp['weight_water'])            

            if not mp.loc[0,'ROMol'].HasSubstructMatch(Chem.MolFromSmarts(smarts)):
                print('Error: cannot find smarts match!\n')
                sys.exit()
            else:
                atom_list=mp.loc[0,'ROMol'].GetSubstructMatches(Chem.MolFromSmarts(smarts))[0]

            df.at[counter,'XTB_Energy']=gather_xtb_ene(weights)
            df.at[counter,'Homo-Lumo']=gather_xtb_hlgap(weights)
            df.at[counter,'Dipole']=gather_xtb_dipole(weights)
            df.at[counter,'EI']=gather_xtb_ei(weights)
            df.at[counter,'EA']=gather_xtb_ea(weights)

            for i in range(0,len(atom_list)):
                f=gather_fukui_ensemble(atom_list[i])
                df.at[counter,'f+_'+str(i)]=Boltzmann_avg(f,weights)

                q=gather_charge_ensemble(atom_list[i],'optimization')
                df.at[counter,'q_'+str(i)]=Boltzmann_avg(q,weights)

                qplus=gather_charge_ensemble(atom_list[i],'qplus')
                df.at[counter,'q+_'+str(i)]=Boltzmann_avg(qplus,weights)

                sasa=gather_sasa_ensemble(atom_list[i])
                df.at[counter,'sasa_'+str(i)]=Boltzmann_avg(sasa,weights)

            os.chdir('../')

            counter+=1

    return df

###############################################################################

#! Assume we execute from the directory of a molpipe run
#! Run XTB and gather descriptors
#! Smarts string must be present in all molecules

###############################################################################

parser = argparse.ArgumentParser(description='Run XTB descriptors')
parser.add_argument('-s','--smarts',dest='s',help='Smarts string of reactive group',type=str,required=True)
args=vars(parser.parse_args())

# Smarts of reactive piece
smarts=args['s']

# Name the output file
output_csv='descriptor_xtb_calc.csv'

# Get all SDF files in directory
all_sdf=glob.glob('*sdf')

# Setup directories
print('Preparing directories ... \n')
setup_dir(all_sdf)

# Run XTB
print('Running XTB ... \n')
run_xtb(all_sdf)

# Gather output
print('Ouput data ... \n')
result_df=gather_xtb(all_sdf,smarts)
result_df.to_csv(output_csv,index=False)

