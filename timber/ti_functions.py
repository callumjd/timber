# timber

import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops
from rdkit.Chem import SDWriter
from .molecule_ff import Molecule_ff,add_rd_bonds,get_rdkit_info
from .geometry import cart_distance
from .align import get_mcs,rms_fit,rigid_coordinate_set
from .ligprep_tools import run_antechamber,Info_Mol2,write_rd_pdb,make_off

##############################################################################

def compare_atom(atm1,atm2,tol=0.1):

    # criteria: position, atomic number, atom type, bond count, hybridization, charge difference < 0.2
    criteria=6

    ctr=0
    if abs(cart_distance(atm2.xyz,atm1.xyz))<tol:
        ctr+=1
    if atm1.atomic_num==atm2.atomic_num:
        ctr+=1
    if atm1.atom_type==atm2.atom_type:
        ctr+=1
    if atm1.bond_count==atm2.bond_count:
        ctr+=1
    if atm1.hybrid==atm2.hybrid:
        ctr+=1
    if abs(atm2.atom_charge-atm1.atom_charge)<0.2:
        ctr+=1

    if ctr==criteria:
        return True
    else:
        return False

def compare_hydrogen_mcs(h_idx1,h_idx2,mol_off1,mol_off2,mcs_atoms1,mcs_atoms2):

    ctr=0
    for at1 in mcs_atoms1:
        if mol_off1.query_bond_atomIdx(h_idx1,at1):
            ctr+=1
            break

    for at2 in mcs_atoms2:
        if mol_off2.query_bond_atomIdx(h_idx2,at2):
            ctr+=1
            break

    if ctr==2:
        return True
    else:
        return False

def compare_mols(rd_mol1,rd_mol2,mol_off1,mol_off2,mcs=None):

    if mcs:
        mcs_atoms1=rd_mol1.GetSubstructMatch(AllChem.MolFromSmarts(mcs))
        mcs_atoms2=rd_mol2.GetSubstructMatch(AllChem.MolFromSmarts(mcs))

    match=[]
    for i in range(0,len(mol_off2.atoms)):
        for j in range(0,len(mol_off1.atoms)):
            if compare_atom(mol_off1.atoms[j],mol_off2.atoms[i]):
                if mcs:
                    if (mol_off2.atoms[i].atomic_num!=1 and mol_off1.atoms[j].atomic_num!=1):
                        if (mol_off2.atoms[i].idx in mcs_atoms2) and (mol_off1.atoms[j].idx in mcs_atoms1):
                            match.append(j)
                    else:
                        # check hydrogen is bonded to mcs atom
                        if compare_hydrogen_mcs(mol_off1.atoms[j].idx,mol_off2.atoms[i].idx,mol_off1,mol_off2,mcs_atoms1,mcs_atoms2):
                            match.append(j)
                else:
                    match.append(j)

    return match

def update_ti_atoms(mol_list,off_list,mcs=None):
    assert len(mol_list)==2
    assert len(off_list)==2

    periodic={'6':'C','1':'H','8':'O','7':'N','17':'Cl','9':'F','16':'S','35':'Br','15':'P','53':'I'}

    for atom in off_list[0].atoms:
        if str(atom.atomic_num) not in periodic:
            raise Exception('Molecule contains element outside of scope')

    for atom in off_list[1].atoms:
        if str(atom.atomic_num) not in periodic:
            raise Exception('Molecule contains element outside of scope')

    matches=compare_mols(mol_list[0],mol_list[1],off_list[0],off_list[1],mcs)

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
                mol_amber.atoms[i].ti_core=True
            elif i in write_last:
                mol_amber.atoms[i].ti_core=False

        for i in write_core:
            new_atom_name=periodic[str(mol_amber.atoms[i].atomic_num)]+str(ele_count[int(mol_amber.atoms[i].atomic_num)])
            mol_amber.atoms[i].name=new_atom_name
            ele_count[int(mol_amber.atoms[i].atomic_num)]+=1

        for i in range(0,len(mol.GetAtoms())):
            if mol_amber.atoms[i].ti_core==False:
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
        if not atom.ti_core:
            ti_region1.append(atom.name)

    ti_str1=''
    for at in ti_region1:
        ti_str1=ti_str1+str(at)+','

    ti_region2=[]
    for atom in off_list[1].atoms:
        if not atom.ti_core:
            ti_region2.append(atom.name)

    ti_str2=''
    for at in ti_region2:
        ti_str2=ti_str2+str(at)+','

    with open(output_file,'w') as f:
        f.write('%s\n' % (ti_str1))
        f.write('%s\n' % (ti_str2))

def run_rbfe_setup(rd_mol1,rd_mol2,dir_1_name='core',dir_2_name='sec_lig',full_mcs=None,align=False):

    if not rd_mol1:
        raise Exception('Error: null mol')
    if not rd_mol2:
        raise Exception('Error: null mol')

    name1=rd_mol1.GetProp('_Name')
    name2=rd_mol2.GetProp('_Name')

    pair_dir=name1+'~'+name2
    if not os.path.exists(pair_dir):
        os.mkdir(pair_dir)
        os.chdir(pair_dir)
    else:
        raise Exception('Error: directory %s exists!' % (pair_dir))

    # Preparing name1 -> name2
    print('%s -> %s \n' % (name1,name2))

    os.mkdir(dir_1_name)
    os.mkdir(dir_2_name)

    # run align if align=True - we expect "full_mcs" to be a smarts string
    local_mcs=get_mcs([Chem.RemoveHs(rd_mol1),Chem.RemoveHs(rd_mol2)],strict=False)
    if align and full_mcs:
        # align lig2 to lig1
        rms_fit(rd_mol1,rd_mol2,mcss=local_mcs.smartsString,mcss_exclusion=full_mcs,bak_seed=full_mcs,tolerance=2.0,ene_cutoff=25)
        rigid_coordinate_set(rd_mol1,rd_mol2,ene_cutoff=35,snap_tol=0.5)

    # write files and create parameters
    parm_mols=[]
    parm_off=[]

    # LIG lig1
    os.chdir(dir_1_name)

    writer=SDWriter('for_parm.sdf')
    writer.write(rd_mol1)
    writer.flush()

    # IMPORTANT - must pass net_charge=Chem.molops.GetFormalCharge(rd_mol1) to get AM1-BCC
    run_antechamber('for_parm.sdf',residue_name='UNL',ff='gaff2')

    LIG=Molecule_ff(name='LIG')
    mol=Chem.SDMolSupplier('UNL.sdf',removeHs=False,sanitize=False)[0]
    mol=AllChem.AssignBondOrdersFromTemplate(rd_mol1,mol) # protect against poor SDF file ...
    Chem.SanitizeMol(mol)
    LIG=get_rdkit_info(mol,LIG)
    LIG=Info_Mol2('UNL.mol2',LIG,len(mol.GetAtoms()),fields=['name','type','charge'])

    parm_mols.append(Chem.Mol(mol))
    parm_off.append(LIG)

    os.chdir('../')

    # MOD lig2
    os.chdir(dir_2_name)

    writer=SDWriter('for_parm.sdf')
    writer.write(rd_mol2)
    writer.flush()

    run_antechamber('for_parm.sdf',residue_name='UNL',ff='gaff2')

    MOD=Molecule_ff(name='MOD')
    mol=Chem.SDMolSupplier('UNL.sdf',removeHs=False,sanitize=False)[0]
    mol=AllChem.AssignBondOrdersFromTemplate(rd_mol2,mol) # protect against poor SDF file ...
    Chem.SanitizeMol(mol)
    MOD=get_rdkit_info(mol,MOD)
    MOD=Info_Mol2('UNL.mol2',MOD,len(mol.GetAtoms()),fields=['name','type','charge'])

    parm_mols.append(Chem.Mol(mol))
    parm_off.append(MOD)

    os.chdir('../')

    # now compare lig1 + lig2, identify soft core atoms
    # pass a new copy of the off objects since they get modified
    # return re-ordered [mol1,mol2] and [off1,off2]
    refit_mols,refit_offs=update_ti_atoms(parm_mols,list(parm_off),mcs=local_mcs.smartsString)

    # LIG refit
    os.chdir(dir_1_name)
    write_rd_pdb(refit_offs[0],refit_mols[0],refit_offs[0].name,'LIG.pdb')
    make_off(refit_offs[0],'make_off.leap')
    os.system('tleap -f make_off.leap>out')
    os.system('rm out make_off.leap')

    writer=SDWriter('LIG.sdf')
    writer.write(refit_mols[0])
    writer.flush()

    os.chdir('../')

    # MOD refit
    os.chdir(dir_2_name)
    write_rd_pdb(refit_offs[1],refit_mols[1],refit_offs[1].name,'MOD.pdb')
    make_off(refit_offs[1],'make_off.leap')
    os.system('tleap -f make_off.leap>out')
    os.system('rm out make_off.leap')

    writer=SDWriter('MOD.sdf')
    writer.write(refit_mols[1])
    writer.flush()

    os.chdir('../')

    # TI masks
    write_ti_strings(refit_offs,'TI_MASKS.dat')

    # pair dir
    os.chdir('../')

