#!/usr/bin/env python
import os
import sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem,SDWriter
import timber

##############################################################################

input_sdf='file.sdf'
ff='gaff2'

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

df=pd.read_csv('mapping.csv')

# prepare all dir lig1->lig2
for index,row in df.iterrows():
    name1=row['Name1']
    name2=row['Name2']
    timber.run_rbfe_setup(all_mols[all_names.index(name1)],all_mols[all_names.index(name2)],ff=ff,full_mcs=full_mcs.smartsString,align=True)

