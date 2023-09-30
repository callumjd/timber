#!/usr/bin/env python
import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdmolops
import timber
from timber.geometry import Coord

filename='lig1.sdf'
residue_name='UNL'
ff='gaff2'

mol=Chem.rdmolfiles.MolFromMolFile(filename,removeHs=False)

# run antechamber
timber.run_antechamber(filename,residue_name=residue_name,ff=ff)##net_charge=int(rdmolops.GetFormalCharge(mol)))

# build mol_ff
parm_mol=timber.Molecule_ff(name=residue_name)

parm_mol=timber.get_rdkit_info(mol,parm_mol)

# save the atom types, charges
parm_mol=timber.Info_Mol2(residue_name+'.mol2',parm_mol,len(parm_mol.atoms),fields=['name','type','charge'])

# write PDB and .off
timber.write_rd_pdb(parm_mol,mol,parm_mol.name,parm_mol.name+'.pdb')

timber.make_off(parm_mol,'make_off.leap')
os.system('tleap -f make_off.leap>out')
os.system('rm make_off.leap out')

# build prmtop
if timber.check_file('missing_'+ff+'.frcmod'):
    timber.build_parm(parm_mol.name,ff,'make_lig.leap',frcmod_file='missing_'+ff+'.frcmod')
else:
    timber.build_parm(parm_mol.name,ff,'make_lig.leap')

os.system('tleap -f make_lig.leap>out')

# add force field terms to the parm_mol
parm_mol=timber.get_parmed_info(parm_mol)

# print energy terms
print(timber.return_energy_decomposition(parm_mol))
print(np.sum(timber.return_energy_decomposition(parm_mol)))

