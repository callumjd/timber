#!/usr/bin/env python
import os
import sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem,SDWriter,rdmolops
import timber

##############################################################################

# input settings
input_sdf='file.sdf'
ff='gaff'
prot='ptp1b_protein.cpptraj.pdb'
protocol='one-step'
ti_repeats=1
#res_seq=int() # abfe runs require amber numbering residue for Boresch restraint

# this is a dict, setting lambda windows
schedule={'complex_ligands':9,  # one-step, three-step vdw, absolute
          'solvent_ligands':9,
          'complex_decharge':5,  # three-step, absolute-three-step
          'solvent_decharge':5,
          'complex_recharge':5,  # three-step
          'solvent_recharge':5,
          'complex_restraint':5} # absolute

# make a lomap mapping file
timber.run_lomap(input_sdf,dir_name='lomap_dir',output_name='mapping')

# load molecules
suppl=Chem.SDMolSupplier(input_sdf,removeHs=False)
all_mols=[m for m in suppl]
all_names=[m.GetProp('_Name') for m in suppl]

suppl_noH=Chem.SDMolSupplier(input_sdf,removeHs=True)
all_mols_noH=[m for m in suppl_noH]

# MCSS for ligand alignment
full_mcs=timber.get_mcs(all_mols_noH,strict=True)

# Align all molecules to first one in the set
for m in all_mols:
    if m.GetProp('_Name')!=all_names[0]:
        timber.rms_fit(all_mols[0],m,mcss=full_mcs.smartsString,mcss_exclusion=None,bak_seed=None,tolerance=2.0,ene_cutoff=500)

# read in the lomap network
df=pd.read_csv('mapping.csv')

# prepare all dir lig1->lig2
for index,row in df.iterrows():
    name1=row['Name1']
    name2=row['Name2']
    timber.run_rbfe_setup(all_mols[all_names.index(name1)],all_mols[all_names.index(name2)],ff=ff,full_mcs=full_mcs.smartsString,align=True)

    # for absolute runs
    #timber.run_abfe_setup(all_mols[all_names.index(name1)],ff=ff)

# get ion numbers for salt concentration. water number is approximated as prot_res*50
prot_chg,prot_res=timber.protein_charge(prot)
pos_ion,neg_ion=timber.return_salt(nwat=prot_res*50,conc=0.15,charge=prot_chg+int(rdmolops.GetFormalCharge(all_mols[0])))

# make build.leap file and build the prmtop
timber.build_ti_leap(prot,ff=ff,protocol=protocol,pos_ion=pos_ion,neg_ion=neg_ion)

if ff=='openff':
    timber.run_build(df,protocol=protocol,hmass=True,use_openff=True)
else:
    timber.run_build(df,protocol=protocol,hmass=True,use_openff=False)

# absolute run: write DISANG file - requires an amber numbered res_seq to set residue index for Boresch restraints
#timber.write_abfe_disang(df,protocol=protocol,res_seq=res_seq)

# submit windows
timber.run_prod(df,protocol=protocol,ti_repeats=ti_repeats,schedule=schedule,hmass=True,equil_ns=1,prod_ns=2,monte_water=1)

# example to extend the prod sampling (using same prod.in file)
timber.extend_prod(df)

# run_analysis once runs are complete
timber.run_analysis(df)

# cinnabar
timber.run_cinnabar(df)

