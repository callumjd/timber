# timber

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign,rdqueries,rdMolTransforms,rdFMCS
from rdkit.Chem import ChemicalForceFields,rdForceFieldHelpers

################################################################################

def compute_msd(m1, m2, matches1, matches2):
    """Compute the lowest MSD between two molecules m1 and m2.
    
    Only atom indices listed in the matches1 and matches2
    vectors of idx vectors are considered.
    
    Return a tuple composed by:
    * best_msd across the molecule pair
    * best_match1: idx vector with atom indices for m1
    * best_match2: idx vector with atom indices for m2
    * best_weights: per-atom weights in the 0.0-1.0 range
      weights, which are inversely proportional to the original distance
      between the corresponding atom pair, can be used to give more weight
      in the subsequent RMS fit to the atom pairs which are closer in
      the original poses
    """
    best_msd = float("inf")
    best_match1 = None
    best_match2 = None
    c1 = m1.GetConformer()
    c2 = m2.GetConformer()
    for match1 in matches1:
        for match2 in matches2:
            n_matches = len(match1)
            assert len(match1) == len(match2)
            sq_distances = np.array([(c1.GetAtomPosition(match1[i]) - c2.GetAtomPosition(match2[i])).LengthSq() for i in range(n_matches)])
            msd = np.average(sq_distances)
            if (msd < best_msd):
                best_msd = msd
                best_match1 = match1
                best_match2 = match2
                best_weights = np.reciprocal(sq_distances)
                min_sq_dist = np.min(best_weights)
                max_sq_dist = np.max(best_weights)
                r = max_sq_dist - min_sq_dist
                if (r < 1.e-5):
                    best_weights.fill(1.0)
                else:
                    best_weights -= min_sq_dist
                    best_weights /= (max_sq_dist - min_sq_dist)
    return best_msd, best_match1, best_match2, best_weights

def get_mcs(mol_list,seed=None,strict=False):

    # strict matching
    if strict:
        if seed:
            mcs=rdFMCS.FindMCS(mol_list,timeout=60,atomCompare=rdFMCS.AtomCompare.CompareElements,bondCompare=rdFMCS.BondCompare.CompareOrder,completeRingsOnly=True,matchValences=True,seedSmarts=seed)
        else:
            mcs=rdFMCS.FindMCS(mol_list,timeout=60,atomCompare=rdFMCS.AtomCompare.CompareElements,bondCompare=rdFMCS.BondCompare.CompareOrder,completeRingsOnly=True,matchValences=True)

    # loose matching
    else:
        if seed:
            mcs=rdFMCS.FindMCS(mol_list,timeout=60,atomCompare=rdFMCS.AtomCompare.CompareAny,bondCompare=rdFMCS.BondCompare.CompareAny,completeRingsOnly=True,matchValences=True,seedSmarts=seed)
        else:
            mcs=rdFMCS.FindMCS(mol_list,timeout=60,atomCompare=rdFMCS.AtomCompare.CompareAny,bondCompare=rdFMCS.BondCompare.CompareAny,completeRingsOnly=True,matchValences=True)

    return mcs

def max_mcss(mol_in,mol_list,seed=None):

    sma_mcss=[]
    len_mcss=[]
    for mol in mol_list:
        res=get_mcs([mol_in,mol],seed=seed,strict=False)
        sma_mcss.append(res)
        len_mcss.append(res.numAtoms)

    if len(len_mcss)>0:
        max_idx=len_mcss.index(max(len_mcss)) 
        return mol_list[max_idx],sma_mcss[max_idx]
    else:
        return mol_list[0],seed

def rigid_coordinate_set(mol1,mol2,ene_cutoff=500,snap_tol=5):

    # zero energy
    mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol2)
    ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol2, mp)
    zero = ff.CalcEnergy()

    if snap_tol<5:
        for r_at in mol1.GetAtoms():
            for f_at in mol2.GetAtoms():
                r_pos=mol1.GetConformer().GetAtomPosition(r_at.GetIdx())
                f_pos=mol2.GetConformer().GetAtomPosition(f_at.GetIdx())
                dist = (f_pos - r_pos).Length()
                if dist<snap_tol:
                    if (r_at.GetAtomicNum()==1) and (f_at.GetAtomicNum()==1):
                        mol2.GetConformer().SetAtomPosition(f_at.GetIdx(),r_pos)
                    elif (r_at.GetAtomicNum()!=1) and (f_at.GetAtomicNum()!=1):
                        mol2.GetConformer().SetAtomPosition(f_at.GetIdx(),r_pos)

                    mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol2)
                    ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol2, mp)
                    ene = ff.CalcEnergy()

                    if (ene-zero)<ene_cutoff:
                        zero = ene
                    else:
                        mol2.GetConformer().SetAtomPosition(f_at.GetIdx(),f_pos)

def rms_fit(mol1,mol2,mcss,mcss_exclusion=None,bak_seed=None,tolerance=1.0,ene_cutoff=500):

    mcss=mcss.replace(':','~').replace(']-',']~').replace('-[','~[').replace(')-',')~').replace('-(','~(').replace('=','~').replace('#0','*')

    matches1 = mol1.GetSubstructMatches(AllChem.MolFromSmarts(mcss),uniquify=False)
    matches2 = mol2.GetSubstructMatches(AllChem.MolFromSmarts(mcss),uniquify=False)

    if (len(matches1)==0 or len(matches2)==0):
        matches1 = mol1.GetSubstructMatches(AllChem.MolFromSmarts(bak_seed),uniquify=False)
        matches2 = mol2.GetSubstructMatches(AllChem.MolFromSmarts(bak_seed),uniquify=False)

    msd, match1, match2, weights = compute_msd(mol1, mol2, matches1, matches2)

    if mcss_exclusion:
        excluded_atoms = mol2.GetSubstructMatches(AllChem.MolFromSmarts(mcss_exclusion))[0]

    # do the MCSS alignment
    frozen_prb_indices = []

    m, wts = (match1,match2),weights
    m = list(reversed(m))
    atom_map = list(zip(*m))
    wts = [min(max(1.e-5, wi), 1.0 - 1.e-5) for wi in wts]

    rdMolAlign.AlignMol(mol2, mol1, atomMap=atom_map, weights=wts)

    # manual set of coords 
    mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol2)
    ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol2, mp)
    zero = ff.CalcEnergy()

    for prb_idx, ref_idx in atom_map:
        ref_pos = mol1.GetConformer().GetAtomPosition(ref_idx)
        prb_pos = mol2.GetConformer().GetAtomPosition(prb_idx)
        dist = (ref_pos - prb_pos).Length()
        if (dist > tolerance):
            pass
        else:
            if mcss_exclusion:
                if prb_idx not in excluded_atoms:
                    mol2.GetConformer().SetAtomPosition(prb_idx, ref_pos)
            else:
                mol2.GetConformer().SetAtomPosition(prb_idx, ref_pos)

        frozen_prb_indices.append(prb_idx)

        # energy check, unless we set ene_cutoff=500
        if ene_cutoff<500:
            mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol2)
            ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol2, mp)
            ene = ff.CalcEnergy()

            if (ene-zero)<ene_cutoff:
                zero = ene
            else:
                mol2.GetConformer().SetAtomPosition(prb_idx, prb_pos)

    # relax atoms outside mcss
    if frozen_prb_indices:
        mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol2)
        ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol2, mp)

        for idx in range(0,len(mol2.GetAtoms())):
            if idx in frozen_prb_indices:
                ff.MMFFAddPositionConstraint(idx, 0.0, 1.e4)
            else:
                ff.MMFFAddPositionConstraint(idx, 0.0, 1.e2)

        ff.Minimize(maxIts=10000, forceTol=0.1)

def coord_min(mol_list,mcss):

    for m in mol_list:
        if m.HasSubstructMatch(AllChem.MolFromSmarts(mcss)):
            freeze=m.GetSubstructMatches(AllChem.MolFromSmarts(mcss))[0]

            mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(m)
            ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(m, mp)

            if ff:
                for i in range(0,len(m.GetAtoms())):
                    if i not in freeze:
                        if int(m.GetAtomWithIdx(i).GetAtomicNum())!=1:
                            ff.MMFFAddPositionConstraint(i, 0.0, 1)
                    else:
                        ff.MMFFAddPositionConstraint(i, 0.0, 1.e6)

                ff.Minimize(maxIts=10000, forceTol=0.1)

