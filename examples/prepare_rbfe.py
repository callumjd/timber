#!/usr/bin/env python
import os
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem,SDWriter
import timber

##############################################################################

input_sdf='ptp1b_ligands.sdf'
name1='lig_20667_2qbp'
name2='lig_20669_2qbr'

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

# example setup for lig1 and lig2
# here, we would iterate over all transform pairs
pair_dir=name1+'~'+name2
os.mkdir(pair_dir)
os.chdir(pair_dir)

lig1=all_mols[all_names.index(name1)]
lig2=all_mols[all_names.index(name2)]

# align lig2 to lig1
local_mcs=timber.get_mcs([all_mols_noH[all_names.index(name1)],all_mols_noH[all_names.index(name2)]],strict=False)

timber.rms_fit(lig1,lig2,mcss=local_mcs.smartsString,mcss_exclusion=full_mcs.smartsString,bak_seed=full_mcs.smartsString,tolerance=2.0,ene_cutoff=25)
timber.rigid_coordinate_set(lig1,lig2,ene_cutoff=35,snap_tol=0.5)

# write files and create parameters
parm_mols=[]
parm_off=[]

# LIG
os.mkdir('core')
os.chdir('core')

writer=SDWriter('for_parm.sdf')
writer.write(lig1)
writer.flush()

timber.run_antechamber('for_parm.sdf',residue_name='UNL',ff='gaff2')

LIG=timber.Molecule_ff(name='LIG')
rd_mol=Chem.SDMolSupplier('UNL.sdf',removeHs=False,sanitize=False)[0]
rd_mol=AllChem.AssignBondOrdersFromTemplate(lig1,rd_mol)
Chem.SanitizeMol(rd_mol)
LIG=timber.get_rdkit_info(rd_mol,LIG)
LIG=timber.Info_Mol2('UNL.mol2',LIG,len(rd_mol.GetAtoms()),fields=['name','type','charge'])

parm_mols.append(Chem.Mol(rd_mol))
parm_off.append(LIG)

os.chdir('../')

# MOD
os.mkdir('sec_lig')
os.chdir('sec_lig')

writer=SDWriter('for_parm.sdf')
writer.write(lig2)
writer.flush()

timber.run_antechamber('for_parm.sdf',residue_name='UNL',ff='gaff2')

MOD=timber.Molecule_ff(name='MOD')
rd_mol=Chem.SDMolSupplier('UNL.sdf',removeHs=False,sanitize=False)[0]
rd_mol=AllChem.AssignBondOrdersFromTemplate(lig2,rd_mol)
Chem.SanitizeMol(rd_mol)
MOD=timber.get_rdkit_info(rd_mol,MOD)
MOD=timber.Info_Mol2('UNL.mol2',MOD,len(rd_mol.GetAtoms()),fields=['name','type','charge'])

parm_mols.append(Chem.Mol(rd_mol))
parm_off.append(MOD)

os.chdir('../')

# now compare lig1 + lig2, identify soft core atoms
# pass a new copy of the off objects since they get modified
# return re-ordered [mol1,mol2] and [off1,off2]
refit_mols,refit_offs=timber.update_ti_atoms(parm_mols,list(parm_off),mcs=local_mcs.smartsString)

# LIG refit
os.chdir('core')
timber.write_rd_pdb(refit_offs[0],refit_mols[0],refit_offs[0].name,'LIG.pdb')
timber.make_off(refit_offs[0],'make_off.leap')
os.system('tleap -f make_off.leap>out')
os.system('rm out make_off.leap')

writer=SDWriter('LIG.sdf')
writer.write(refit_mols[0])
writer.flush()

os.chdir('../')

# MOD refit
os.chdir('sec_lig')
timber.write_rd_pdb(refit_offs[1],refit_mols[1],refit_offs[1].name,'MOD.pdb')
timber.make_off(refit_offs[1],'make_off.leap')
os.system('tleap -f make_off.leap>out')
os.system('rm out make_off.leap')

writer=SDWriter('MOD.sdf')
writer.write(refit_mols[1])
writer.flush()

os.chdir('../')

# TI masks
timber.write_ti_strings(refit_offs,'TI_MASKS.dat')

# pair dir
os.chdir('../')

