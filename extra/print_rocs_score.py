#!/usr/bin/env python
from openeye import oechem
from openeye import oeshape
from openeye.oechem import *

ref_file='ref.sdf'
sdf_file='probe.sdf'

# Read mol
ifs=oemolistream()

refmol=OEMol()
ifs=oemolistream(ref_file)
OEReadMolecule(ifs, refmol)

sdf_fs=oemolistream()
if sdf_fs.open(sdf_file):
    for fitmol in sdf_fs.GetOEGraphMols():

        prep=oeshape.OEOverlapPrep()
        prep.Prep(refmol)
        func=oeshape.OEOverlapFunc()
        func.SetupRef(refmol)

        res=oeshape.OEOverlapResults()

        prep.Prep(fitmol)
        func.Overlap(fitmol, res)

        print('Combo: %s %lf' % (fitmol.GetTitle(),float(res.GetTanimotoCombo())))
        print('Shape: %s %lf' % (fitmol.GetTitle(),float(res.GetTanimoto())))
        print('Colour: %s %lf' % (fitmol.GetTitle(),float(res.GetColorTanimoto())))

