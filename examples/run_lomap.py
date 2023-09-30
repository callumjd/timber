#!/usr/bin/env python
import os
import lomap
from rdkit import Chem
from rdkit.Chem import SDWriter

input_sdf='ptp1b_ligands.sdf'

# load molecules
suppl=Chem.SDMolSupplier(input_sdf,removeHs=False)

if not os.path.isdir('lomap_dir'):
    os.mkdir('lomap_dir')

for m in suppl:
    name=m.GetProp('_Name')
    writer=SDWriter('lomap_dir/'+name+'.sdf')
    writer.write(m)
    writer.flush()

# allow charge changing with ecrscore
db_mol = lomap.DBMolecules('lomap_dir/', output=True, ecrscore=0)

strict, loose = db_mol.build_matrices()

nx_graph = db_mol.build_graph()


