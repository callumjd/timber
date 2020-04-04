#!/usr/bin/env python

#
# TIMBER
# Callum J Dickson
# 7 March 2020
#

import argparse
import os
import sys
import re
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import SDWriter
from rdkit.Chem import ChemicalForceFields
from rdkit.Chem import BRICS
from rdkit.Chem import rdmolops
import distutils.spawn
from ligprep_tools import *

###############################################################################

def check_file(file_in):
    output=False
    if os.path.exists(file_in) and os.path.getsize(file_in)>0:
        output=True
    return output

def mapping_tuples(csv_file):
    map_list=[]

    # assumes first line are comments
    with open(csv_file,'r') as f:
        f.readline()
        for line in f:
            lig1=line.split(',')[0].strip()
            lig2=line.split(',')[1].strip()
            map_list.append((lig1,lig2))
            
    return map_list

def update_atom_position(mol1,mol2):
    mol_copy=Chem.Mol(mol2)

    # This is a work-around to get a seedSmarts for the FMCS algorithm
    # and prevent the occassional hanging of FMCS
    # Might be unnecessary with future versions of rdkit
    core_frags=BRICS.BRICSDecompose(Chem.RemoveHs(mol1))
    frag_smarts=[]
    for frag in enumerate(core_frags):
        smi_str=(re.sub('[[1-9][0-9]{0,2}\*]','[*]',frag[1]))
        frag_smarts.append(Chem.MolToSmarts(Chem.MolFromSmiles(smi_str)).replace(':','~').replace('-','~').replace('=','~').replace('#0','*'))

    seed=None
    seed_hits=[]
    for query in frag_smarts:
        if mol_copy.HasSubstructMatch(Chem.MolFromSmarts(query)):
            seed_hits.append(query)

    seed_hits.sort(key=len,reverse=True)
    seed=seed_hits[0]

    # Now get MCSS
    res=rdFMCS.FindMCS([mol1,mol_copy],seedSmarts=seed)
    mcs_q=Chem.MolFromSmarts(res.smartsString)

    # Get atom IDs
    template_list=mol1.GetSubstructMatches(mcs_q)
    hit_atom_list=mol_copy.GetSubstructMatches(mcs_q)

    ss_match_pairs=[]
    ss_distance=[]
    for template_opts in template_list:
        for hit_atom_opts in hit_atom_list:
            ss_match_pairs.append((template_opts,hit_atom_opts))
            running_distance=0
            for i in range(len(template_opts)):
                origin=mol1.GetConformer().GetAtomPosition(template_opts[i])
                pos=mol_copy.GetConformer().GetAtomPosition(hit_atom_opts[i])

                p1=np.array([origin.x,origin.y,origin.z])
                p2=np.array([pos.x,pos.y,pos.z])

                sq_dist=np.sum((p1-p2)**2,axis=0)
                dist=np.sqrt(sq_dist)

                running_distance+=dist
            ss_distance.append(running_distance)

    index_min_pair=ss_distance.index(min(ss_distance))
    template,hit_atom=ss_match_pairs[index_min_pair]

    # Update XYZ coords of MCSS
    running_distance=0
    for i in range(0,len(template)):
        origin=mol1.GetConformer().GetAtomPosition(template[i])
        pos=mol_copy.GetConformer().GetAtomPosition(hit_atom[i])

        p1=np.array([origin.x,origin.y,origin.z])
        p2=np.array([pos.x,pos.y,pos.z])

        sq_dist=np.sum((p1-p2)**2,axis=0) 
        dist=np.sqrt(sq_dist)

        running_distance+=dist

        if mol_copy.GetAtomWithIdx(hit_atom[i]).GetAtomicNum()!=1:
            mol_copy.GetConformer().SetAtomPosition(hit_atom[i],(origin.x,origin.y,origin.z))

    if running_distance>0.1:
        # relax atoms outside MCSS
        res_atom=[]
        for atom in mol_copy.GetAtoms():
            if atom.GetIdx() not in hit_atom:
                res_atom.append(atom.GetIdx())

        # do minimization
        mp=ChemicalForceFields.MMFFGetMoleculeProperties(mol_copy)
        ff=ChemicalForceFields.MMFFGetMoleculeForceField(mol_copy,mp)

        for val in hit_atom:
            if mol_copy.GetAtomWithIdx(val).GetAtomicNum()!=1:
                ff.AddFixedPoint(val)
        for val in res_atom:
            ff.MMFFAddPositionConstraint(val,1,5)

        ff.Minimize()

    # repeat again to fix hydrogen positions
    ss_match_pairs=[]
    ss_distance=[]
    for template_opts in template_list:
        for hit_atom_opts in hit_atom_list:
            ss_match_pairs.append((template_opts,hit_atom_opts))
            running_distance=0
            for i in range(len(template_opts)):
                origin=mol1.GetConformer().GetAtomPosition(template_opts[i])
                pos=mol_copy.GetConformer().GetAtomPosition(hit_atom_opts[i])

                p1=np.array([origin.x,origin.y,origin.z])
                p2=np.array([pos.x,pos.y,pos.z])

                sq_dist=np.sum((p1-p2)**2,axis=0)
                dist=np.sqrt(sq_dist)

                running_distance+=dist
            ss_distance.append(running_distance)

    index_min_pair=ss_distance.index(min(ss_distance))
    template,hit_atom=ss_match_pairs[index_min_pair]

    # Update XYZ coords of MCSS
    running_distance=0
    for i in range(0,len(template)):
        origin=mol1.GetConformer().GetAtomPosition(template[i])
        pos=mol_copy.GetConformer().GetAtomPosition(hit_atom[i])

        p1=np.array([origin.x,origin.y,origin.z])
        p2=np.array([pos.x,pos.y,pos.z])

        sq_dist=np.sum((p1-p2)**2,axis=0)
        dist=np.sqrt(sq_dist)

        running_distance+=dist

        mol_copy.GetConformer().SetAtomPosition(hit_atom[i],(origin.x,origin.y,origin.z))

    return mol_copy

## May want to relax atom_type comparison here using GAFF or GAFF2 ##
def compare_atom(atm1,atm2,tol=0.1):
    if (atm1.atomic_num==atm2.atomic_num) and (atm1.atom_type==atm2.atom_type) and (atm1.hybrid==atm2.hybrid) and (atm1.bond_count==atm2.bond_count) and abs(atm1.x-atm2.x)<tol and abs(atm1.y-atm2.y)<tol and abs(atm1.z-atm2.z)<tol:
        return True
    else:
        return False

def compare_mols(mol_off1,mol_off2):
    match=[]
    for i in range(0,len(mol_off2.atoms)):
        for j in range(0,len(mol_off1.atoms)):
            if compare_atom(mol_off1.atoms[j],mol_off2.atoms[i]):
                match.append(j)
    return match

def update_ti_atoms(mol_list,off_list):
    assert len(mol_list)==2
    assert len(off_list)==2

    periodic={'6':'C','1':'H','8':'O','7':'N','17':'Cl','9':'F','16':'S','35':'Br','15':'P','53':'I'}

    matches=compare_mols(off_list[0],off_list[1])

    MCS_atoms_amber=Molecule_ff(name='MCS_atoms')
    for i in matches:
        MCS_atoms_amber.add_atom(off_list[0].atoms[i])

    out_mols=[]
    out_off=[]
    for mol,mol_amber in zip(mol_list,off_list):
        ele_count=dict([(6,1),(1,1),(8,1),(7,1),(17,1),(9,1),(16,1),(35,1),(15,1),(53,1)])

        write_core=[]
        write_last=[]

        mol_copy=Chem.Mol(mol)

        for i in range(0,len(MCS_atoms_amber.atoms)):
            for j in range(0,len(mol.GetAtoms())):
                if compare_atom(MCS_atoms_amber.atoms[i],mol_amber.atoms[j]) and j not in write_core:
                    write_core.append(j)

        for i in range(0,len(mol.GetAtoms())):
            if i not in write_core:
                write_last.append(i)

        for i in range(0,len(mol.GetAtoms())):
            if i in write_core:
                mol_amber.atoms[i].core=True
            elif i in write_last:
                mol_amber.atoms[i].core=False

        for i in write_core:
            new_atom_name=periodic[str(mol_amber.atoms[i].atomic_num)]+str(ele_count[int(mol_amber.atoms[i].atomic_num)])
            mol_amber.atoms[i].name=new_atom_name
            ele_count[int(mol_amber.atoms[i].atomic_num)]+=1

        for i in range(0,len(mol.GetAtoms())):
            if mol_amber.atoms[i].core==False:
                new_atom_name=periodic[str(mol_amber.atoms[i].atomic_num)]+str(ele_count[int(mol_amber.atoms[i].atomic_num)])
                mol_amber.atoms[i].name=new_atom_name
                ele_count[int(mol_amber.atoms[i].atomic_num)]+=1

        # return a re-ordered mol
        mol_copy=rdmolops.RenumberAtoms(mol_copy,write_core+write_last)
        out_mols.append(mol_copy)

        # return matching re-ordered amber off
        mol_amber.reorder_atoms(write_core+write_last)
        mol_amber.clear_bonds()
        mol_amber=add_rd_bonds(mol_copy,mol_amber)
        out_off.append(mol_amber)

    return out_mols,out_off

def write_ti_strings(off_list,output_file):
    ti_region1=[]
    for atom in off_list[0].atoms:
        if not atom.core:
            ti_region1.append(atom.name)

    ti_str1=''
    for at in ti_region1:
        ti_str1=ti_str1+str(at)+','

    ti_region2=[]
    for atom in off_list[1].atoms:
        if not atom.core:
            ti_region2.append(atom.name)

    ti_str2=''
    for at in ti_region2:
        ti_str2=ti_str2+str(at)+','

    with open(output_file,'w') as f:
        f.write('%s\n' % (ti_str1))
        f.write('%s\n' % (ti_str2))

###############################################################################

if __name__=='__main__':
    ## Command line arguments ##
    parser = argparse.ArgumentParser(description='TIMBER code for AMBER TI setup\n')
    parser.add_argument('-i',help='CSV file with ligand mappings',required=True)
    parser.add_argument('-sdf',help='SDF file with ligands',required=False)
    parser.add_argument('-ff',help='Force field',choices=['gaff','gaff2'],required=False)
    parser.add_argument('-m',help='Mode',choices=['setup'],required=True)

    args=vars(parser.parse_args())

###############################################################################
## MODE: SETUP ##
###############################################################################

    ## Check files ##
    if args['m']=='setup':
        if args['i']!=None and check_file(args['i']):
            map_file=args['i']
        else:
            print('Error: CSV file required.\n')
            sys.exit()

        if args['sdf']!=None and check_file(args['sdf']):
            ligand_sdf=args['sdf']
        else:
            print('Error: setup mode requires input SDF ligand file.\n')
            sys.exit()

        if args['ff']!=None:
            ff=args['ff']
        else:
            print('Error: setup mode required force field specification.\n')
            sys.exit()

        if not distutils.spawn.find_executable('antechamber'):
            print('Error: cannot find antechamer.\n')
            sys.exit()

    ## Proceed with setup ##
        print('\nSetup: writing transformation directories.\n')
   
        dir_1_name='start'
        dir_2_name='endpoint'

    ## Get a list of the mapping tuples ##
        map_list=mapping_tuples(map_file)

    ## Load the ligand file and save names ##
        ligands=Chem.SDMolSupplier(ligand_sdf,removeHs=False)
        ligands_name=[]
        for mol in ligands:
            ligands_name.append(mol.GetProp('_Name'))

    ## Make directories for each transformation ##
        for pair in map_list:
            print('%s -> %s \n' % (pair[0],pair[1]))

            pair_dir=pair[0]+'_'+pair[1]
            if not os.path.isdir(pair_dir):
                os.mkdir(pair_dir)
                os.chdir(pair_dir)

                os.mkdir(dir_1_name)
                os.mkdir(dir_2_name)

    ## Write start ligand file, parameters ##
                os.chdir(dir_1_name)
                writer=SDWriter('for_parm.sdf')
                if ligands_name.count(pair[0])>0:
                    writer.write(ligands[ligands_name.index(pair[0])])
                    writer.flush()
                else:
                    print('Error: cannot map ligand %s.\n' % (pair[0]))
                    sys.exit()
                run_antechamber('for_parm.sdf','UNL',ff,int(rdmolops.GetFormalCharge(ligands[ligands_name.index(pair[0])])),clean_sdf=True)

                # setup Molecule_ff
                LIG=Molecule_ff(name='LIG')
                n_atoms=len(ligands[ligands_name.index(pair[0])].GetAtoms())
                for at in ligands[ligands_name.index(pair[0])].GetAtoms():
                    x=ligands[ligands_name.index(pair[0])].GetConformer().GetAtomPosition(at.GetIdx()).x
                    y=ligands[ligands_name.index(pair[0])].GetConformer().GetAtomPosition(at.GetIdx()).y
                    z=ligands[ligands_name.index(pair[0])].GetConformer().GetAtomPosition(at.GetIdx()).z
                    LIG.add_atom(Atom_ff(idx=at.GetIdx(),atomic_num=at.GetAtomicNum(),atomic_weight=Chem.GetPeriodicTable().GetAtomicWeight(at.GetAtomicNum()),hybrid=at.GetHybridization(),bond_count=len(at.GetBonds()),x=x,y=y,z=z))
                LIG=Info_Mol2('UNL.mol2',LIG,n_atoms,fields=['name','type','charge'])
                LIG=add_rd_bonds(ligands[ligands_name.index(pair[0])],LIG)

                os.chdir('../')

    ## Write endpoint ligand file, parameters ##
                os.chdir(dir_2_name)
                writer=SDWriter('for_parm.sdf')
                if ligands_name.count(pair[1])>0:
    ## Check and fix XYZ coords of transform ligand ##
                    fix_mol=update_atom_position(ligands[ligands_name.index(pair[0])],ligands[ligands_name.index(pair[1])])

                    writer.write(fix_mol)
                    writer.flush()
                else:
                    print('Error: cannot map ligand %s.\n' % (pair[1]))
                    sys.exit()
                run_antechamber('for_parm.sdf','UNL',ff,int(rdmolops.GetFormalCharge(ligands[ligands_name.index(pair[1])])),clean_sdf=True)

                # setup Molecule_ff
                MOD=Molecule_ff(name='MOD')
                n_atoms=len(fix_mol.GetAtoms())
                for at in fix_mol.GetAtoms(): 
                    x=fix_mol.GetConformer().GetAtomPosition(at.GetIdx()).x
                    y=fix_mol.GetConformer().GetAtomPosition(at.GetIdx()).y
                    z=fix_mol.GetConformer().GetAtomPosition(at.GetIdx()).z
                    MOD.add_atom(Atom_ff(idx=at.GetIdx(),atomic_num=at.GetAtomicNum(),atomic_weight=Chem.GetPeriodicTable().GetAtomicWeight(at.GetAtomicNum()),hybrid=at.GetHybridization(),bond_count=len(at.GetBonds()),x=x,y=y,z=z))
                MOD=Info_Mol2('UNL.mol2',MOD,n_atoms,fields=['name','type','charge'])
                MOD=add_rd_bonds(fix_mol,MOD)

                os.chdir('../')

    ## Now rename and re-order start and endpoint ligand atoms so that TI region is at the end
                parm_mols=[]
                parm_off=[LIG,MOD]
                for parm_dir in [dir_1_name,dir_2_name]:
                    mol=Chem.SDMolSupplier(parm_dir+'/UNL.sdf',removeHs=False,sanitize=False)[0]
                    parm_mols.append(mol)

                # pass a new copy of the off objects since they get modified
                # return re-ordered [mol1,mol2] and [off1,off2]
                refit_mols,refit_offs=update_ti_atoms(parm_mols,list(parm_off))

                os.chdir(dir_1_name)
                write_rd_pdb(refit_offs[0],refit_mols[0],refit_offs[0].name,'LIG.pdb')
                make_off(refit_offs[0],'make_off.leap')
                os.system('tleap -f make_off.leap>out')
                os.system('rm out make_off.leap')
                os.chdir('../')

                os.chdir(dir_2_name)
                write_rd_pdb(refit_offs[1],refit_mols[1],refit_offs[1].name,'MOD.pdb')
                make_off(refit_offs[1],'make_off.leap')
                os.system('tleap -f make_off.leap>out')
                os.system('rm out make_off.leap')
                os.chdir('../')
        
                write_ti_strings(refit_offs,'TI_MASKS.dat')

    ## Exit pair directory
                os.chdir('../')

            else:
                print('Error: directory exists %s.\n' % (pair_dir))

    print('Setup complete.\n')

################################################################################
