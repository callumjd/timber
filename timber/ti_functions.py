# timber

import os
import glob as glob
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops
from rdkit.Chem import SDWriter
from .molecule_ff import Molecule_ff,add_rd_bonds,get_rdkit_info
from .geometry import cart_distance
from .align import get_mcs,rms_fit,rigid_coordinate_set
from .ligprep_tools import run_antechamber,Info_Mol2,write_rd_pdb,make_off,check_file,setup_hmass
from .ti_inputs import write_ti_cluster_script,write_ti_inputs
from .openff_ligprep import off_prmtop_converter 

##############################################################################

# this is a dict, setting lambda windows which can be passed to run_prod 
#
# schedule={'complex_ligands':9,  # one-step, three-step vdw, absolute
#          'solvent_ligands':9,
#          'complex_decharge':5,  # three-step, absolute-three-step
#          'solvent_decharge':5,
#          'complex_recharge':5,  # three-step
#          'solvent_recharge':5,
#          'complex_restraint':5} # absolute

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

def run_abfe_setup(rd_mol1,ff='gaff2',dir_1_name='core'):

    if not rd_mol1:
        raise Exception('Error: null mol')

    name1=rd_mol1.GetProp('_Name')

    pair_dir=name1
    if not os.path.exists(pair_dir):
        os.mkdir(pair_dir)
        os.chdir(pair_dir)
    else:
        raise Exception('Error: directory %s exists!' % (pair_dir))

    # Preparing name1 -> name2
    print('%s -> Nothing \n' % (name1))

    os.mkdir(dir_1_name)
    os.chdir(dir_1_name)

    writer=SDWriter('for_parm.sdf')
    writer.write(rd_mol1)
    writer.flush()

    # IMPORTANT - must pass net_charge=Chem.rdmolops.GetFormalCharge(rd_mol1) to get AM1-BCC
    run_antechamber('for_parm.sdf',residue_name='LIG',ff=ff)

    LIG=Molecule_ff(name='LIG')
    mol=Chem.SDMolSupplier('LIG.sdf',removeHs=False,sanitize=False)[0]
    mol=AllChem.AssignBondOrdersFromTemplate(rd_mol1,mol) # protect against poor SDF file ...
    Chem.SanitizeMol(mol)
    LIG=get_rdkit_info(mol,LIG)
    LIG=Info_Mol2('LIG.mol2',LIG,len(mol.GetAtoms()),fields=['name','type','charge'])
    write_rd_pdb(LIG,mol,'LIG','LIG.pdb')
    make_off(LIG,'make_off.leap')
    os.system('tleap -f make_off.leap>out')
    os.system('rm out make_off.leap')

    os.chdir('../../')

def run_rbfe_setup(rd_mol1,rd_mol2,ff='gaff2',dir_1_name='core',dir_2_name='sec_lig',full_mcs=None,align=False):

    if ff=='openff':
        ff='gaff'

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
    # ene_cutoff and snap_tol may need to be adjusted - expose?
    local_mcs=get_mcs([Chem.RemoveHs(rd_mol1),Chem.RemoveHs(rd_mol2)],strict=False)
    if align and full_mcs:
        # align lig2 to lig1
        rms_fit(rd_mol1,rd_mol2,mcss=local_mcs.smartsString,mcss_exclusion=full_mcs,bak_seed=full_mcs,tolerance=2.0,ene_cutoff=75)
        rigid_coordinate_set(rd_mol1,rd_mol2,ene_cutoff=75,snap_tol=0.75)

    # write files and create parameters
    parm_mols=[]
    parm_off=[]

    # LIG lig1
    os.chdir(dir_1_name)

    writer=SDWriter('for_parm.sdf')
    writer.write(rd_mol1)
    writer.flush()

    # IMPORTANT - must pass net_charge=Chem.rdmolops.GetFormalCharge(rd_mol1) to get AM1-BCC
    run_antechamber('for_parm.sdf',residue_name='UNL',ff=ff)

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

    run_antechamber('for_parm.sdf',residue_name='UNL',ff=ff)

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

def check_leap_build(leap_output_file):
    # TO DO: Update with better error checks
    output=True
    with open(leap_output_file,'r') as f:
        data=f.readlines()

    for line in data:
        if 'Parameter file was not saved.' in line:
            output=False
            break
        elif 'Fatal Error' in line:
            output=False
            break
        elif 'Exiting LEaP: Errors' in line:
            err_count=int(line.split()[4].strip(';'))
            if err_count>0:
                output=False
                break

    return output

def run_build(df,protocol='one-step',hmass=True,use_openff=False,dir_1_name='core',dir_2_name='sec_lig'):

    # check if build.leap exists
    if not check_file('build.leap'):
        raise Exception('Error: build mode requires build.leap file\n')

    # build all dir lig1->lig2
    for index,row in df.iterrows():
        name1=row['Name1']
        # rbfe
        if 'Name2' in list(df.columns):
            name2=row['Name2']
            pair_dir=name1+'~'+name2
        # absolute
        else:
            pair_dir=name1

        os.chdir(pair_dir)
        os.mkdir('complex')
        os.mkdir('solvent')

        os.system('cp ../build.leap .')
        os.system('tleap -f build.leap>leap_out')

        if (check_file('complex/complex_ligands.prmtop') and check_file('solvent/solvent_ligands.prmtop') and check_leap_build('leap_out')):
            print('System built: %s\n' % (pair_dir))
            os.system('rm leap_out')
        else:
            raise Exception('Error: could not build system %s\n' % (pair_dir))

        # first, do the openff setup
        if use_openff:
            for media in ['complex','solvent']:
                os.chdir(media)

                if protocol=='one-step':
                    off_prmtop_converter(media+'_ligands',['../'+dir_1_name+'/LIG','../'+dir_2_name+'/MOD'],xml='openff_unconstrained-2.0.0.offxml')
                elif protocol=='three-step':
                    off_prmtop_converter(media+'_ligands',['../'+dir_1_name+'/LIG','../'+dir_2_name+'/MOD'],xml='openff_unconstrained-2.0.0.offxml')
                    off_prmtop_converter(media+'_decharge',['../'+dir_1_name+'/LIG','../'+dir_1_name+'/LIG'],xml='openff_unconstrained-2.0.0.offxml')
                    off_prmtop_converter(media+'_recharge',['../'+dir_2_name+'/MOD','../'+dir_2_name+'/MOD'],xml='openff_unconstrained-2.0.0.offxml')
                elif protocol=='absolute':
                    off_prmtop_converter(media+'_ligands',['../'+dir_1_name+'/LIG'],xml='openff_unconstrained-2.0.0.offxml')
                    if media=='complex':
                        off_prmtop_converter(media+'_restraint',['../'+dir_1_name+'/LIG','../'+dir_1_name+'/LIG'],xml='openff_unconstrained-2.0.0.offxml')

                elif protocol=='absolute-three-step':
                    off_prmtop_converter(media+'_ligands',['../'+dir_1_name+'/LIG'],xml='openff_unconstrained-2.0.0.offxml')
                    off_prmtop_converter(media+'_decharge',['../'+dir_1_name+'/LIG','../'+dir_1_name+'/LIG'],xml='openff_unconstrained-2.0.0.offxml')
                    if media=='complex':
                        off_prmtop_converter(media+'_restraint',['../'+dir_1_name+'/LIG','../'+dir_1_name+'/LIG'],xml='openff_unconstrained-2.0.0.offxml')

                os.chdir('../')

        # prepare the hmass prmtop
        for media in ['complex','solvent']:
            os.chdir(media)
            for prmtop in glob.glob('*prmtop'):
                setup_hmass(prmtop)
            os.chdir('../')

        os.chdir('../')

def build_ti_leap(prot,prep_files=None,pdb_files=None,protocol='one-step',pos_ion=0,neg_ion=0,ff='gaff2',protein_ff='ff19SB',water_ff='tip3p',ion_ff='ionsjc_tip3p',dir_1_name='core',dir_2_name='sec_lig'):

    convert_prep={'off':'loadoff','lib':'loadoff','prep':'loadamberprep','frcmod':'loadamberparams','mol2':'loadmol2','zinc':'source','add':'loadamberparams'}

    if not check_file(prot):
        raise Exception('Error: cannot find %s\n' % (prot))

    with open('build.leap','w') as f:
        if ff=='openff':
            f.write('source leaprc.protein.ff14SB\n')
        else:
            f.write('source leaprc.protein.%s\n' % (protein_ff))
        f.write('source leaprc.water.%s\n' % (water_ff))
        f.write('source leaprc.phosaa19SB\n')
        # DNA force field
        f.write('source leaprc.DNA.bsc1\n')
        f.write('loadamberparams frcmod.%s\n' % (ion_ff))
        # this is for MG ions etc
        f.write('loadamberparams frcmod.ions234lm_126_tip3p\n')

        if ff=='gaff' or ff=='openff':
            f.write('source leaprc.gaff\n')
            f.write('\n')
            f.write('loadamberparams %s/missing_gaff.frcmod\n' % (dir_1_name))
            if protocol not in ['absolute','absolute-three-step']:
                f.write('loadamberparams %s/missing_gaff.frcmod\n' % (dir_2_name))
        elif ff=='gaff2':
            f.write('source leaprc.gaff2\n')
            f.write('\n')
            f.write('loadamberparams %s/missing_gaff2.frcmod\n' % (dir_1_name))
            if protocol not in ['absolute','absolute-three-step']:
                f.write('loadamberparams %s/missing_gaff2.frcmod\n' % (dir_2_name))

        if prep_files and len(prep_files)>0:
            for val in prep_files:
                my_type=val.split('.')[-1]
                f.write('%s %s\n' % (convert_prep[my_type],os.path.abspath(val)))

        f.write('\n')
        f.write('loadoff %s/LIG.off\n' % (dir_1_name))
        if protocol not in ['absolute','absolute-three-step']:
            f.write('loadoff %s/MOD.off\n' % (dir_2_name))
        f.write('\n')
        f.write('LIG=loadpdb %s/LIG.pdb\n' % (dir_1_name))
        if protocol not in ['absolute','absolute-three-step']:
            f.write('MOD=loadpdb %s/MOD.pdb\n' % (dir_2_name))
        f.write('\n')
        f.write('prot=loadpdb %s\n' % os.path.abspath(prot))

        pdb_str=''
        not_wat=''
        if pdb_files:
            for i in range(0,len(pdb_files)):
                f.write('mol%d=loadpdb %s\n' % (i,os.path.abspath(pdb_files[i])))
                pdb_str=pdb_str+'mol'+str(i)+' '
                if not is_water(pdb_files[i]):
                    not_wat=not_wat+'mol'+str(i)+' '

        f.write('\n')

        # initial units
        if protocol in ['one-step','three-step']:
            f.write('ligands = combine {LIG LIG MOD MOD}\n')
        elif protocol in ['absolute','absolute-three-step']:
            f.write('ligands = combine {LIG LIG}\n')

        f.write('solvatebox ligands TIP3PBOX 20.0 0.75\n')
        f.write('ligands_all=copy ligands\n')
        f.write('\n')

        if protocol in ['one-step','three-step']:
            f.write('complex=combine {LIG LIG MOD MOD prot %s}\n' % (pdb_str))
        elif protocol in ['absolute','absolute-three-step']:
            f.write('complex=combine {LIG LIG prot %s}\n' % (pdb_str))

        f.write('solvatebox complex TIP3PBOX 12.0 0.75\n')
        f.write('complex_all=copy complex\n')
        f.write('\n')

        # LIGAND
        # save base ligand prmtops
        if protocol in ['one-step','three-step']:
            f.write('remove ligands ligands.4\n')
            f.write('remove ligands ligands.2\n')
        elif protocol in ['absolute','absolute-three-step']:
            f.write('remove ligands ligands.2\n')
        f.write('addIonsRand ligands Na+ 0\n')
        f.write('addIonsRand ligands Cl- 0\n')
        f.write('savepdb ligands solvent/solvent_ligands.pdb\n')
        f.write('saveamberparm ligands solvent/solvent_ligands.prmtop solvent/solvent_ligands.inpcrd\n')

        if protocol=='three-step':
            f.write('\n')
            f.write('ligands=copy ligands_all\n')
            f.write('remove ligands ligands.4\n')
            f.write('remove ligands ligands.3\n')
            f.write('addIonsRand ligands Na+ 0\n')
            f.write('addIonsRand ligands Cl- 0\n')
            f.write('savepdb ligands solvent/solvent_decharge.pdb\n')
            f.write('saveamberparm ligands solvent/solvent_decharge.prmtop solvent/solvent_decharge.inpcrd\n')
            f.write('\n')
            f.write('ligands=copy ligands_all\n')
            f.write('remove ligands ligands.2\n')
            f.write('remove ligands ligands.1\n')
            f.write('addIonsRand ligands Na+ 0\n')
            f.write('addIonsRand ligands Cl- 0\n')
            f.write('savepdb ligands solvent/solvent_recharge.pdb\n')
            f.write('saveamberparm ligands solvent/solvent_recharge.prmtop solvent/solvent_recharge.inpcrd\n')

        # absolute solvent does not require restraint leg
        elif protocol=='absolute-three-step':
            f.write('\n')
            f.write('ligands=copy ligands_all\n')
            f.write('addIonsRand ligands Na+ 0\n')
            f.write('addIonsRand ligands Cl- 0\n')
            f.write('savepdb ligands solvent/solvent_decharge.pdb\n')
            f.write('saveamberparm ligands solvent/solvent_decharge.prmtop solvent/solvent_decharge.inpcrd\n')

        f.write('\n')

        # COMPLEX
        # save base ligand prmtops
        if protocol in ['one-step','three-step']:
            f.write('remove complex complex.4\n')
            f.write('remove complex complex.2\n')
        elif protocol in ['absolute','absolute-three-step']:
            f.write('remove complex complex.2\n')
        if pos_ion>0 or neg_ion>0:
            f.write('addIonsRand complex Na+ %d\n' % (pos_ion))
            f.write('addIonsRand complex Cl- %d\n' % (neg_ion))
        else:
            f.write('addIonsRand complex Na+ 0\n')
            f.write('addIonsRand complex Cl- 0\n')
        f.write('savepdb complex complex/complex_ligands.pdb\n')
        f.write('saveamberparm complex complex/complex_ligands.prmtop complex/complex_ligands.inpcrd\n')

        if protocol=='three-step':
            f.write('\n')
            f.write('complex=copy complex_all\n')
            f.write('remove complex complex.4\n')
            f.write('remove complex complex.3\n')
            if pos_ion==0 and neg_ion==0:
                f.write('addIonsRand complex Na+ 0\n')
                f.write('addIonsRand complex Cl- 0\n')
            f.write('savepdb complex complex/complex_decharge.pdb\n')
            f.write('saveamberparm complex complex/complex_decharge.prmtop complex/complex_decharge.inpcrd\n')
            f.write('\n')
            f.write('complex=copy complex_all\n')
            f.write('remove complex complex.2\n')
            f.write('remove complex complex.1\n')
            if pos_ion==0 and neg_ion==0:
                f.write('addIonsRand complex Na+ 0\n')
                f.write('addIonsRand complex Cl- 0\n')
            f.write('savepdb complex complex/complex_recharge.pdb\n')
            f.write('saveamberparm complex complex/complex_recharge.prmtop complex/complex_recharge.inpcrd\n')

        elif protocol in ['absolute','absolute-three-step']:
            f.write('\n')
            f.write('complex=copy complex_all\n')
            if pos_ion>0 or neg_ion>0:
                f.write('addIonsRand complex Na+ %d\n' % (pos_ion))
                f.write('addIonsRand complex Cl- %d\n' % (neg_ions))
            else:
                f.write('addIonsRand complex Na+ 0\n')
                f.write('addIonsRand complex Cl- 0\n')
            f.write('savepdb complex complex/complex_restraint.pdb\n')
            f.write('saveamberparm complex complex/complex_restraint.prmtop complex/complex_restraint.inpcrd\n')
            if protocol=='absolute-three-step':
                f.write('savepdb complex complex/complex_decharge.pdb\n')
                f.write('saveamberparm complex complex/complex_decharge.prmtop complex/complex_decharge.inpcrd\n')

        f.write('\n')
        f.write('quit\n')

def set_lambda_values(lambda_input):
    output_lambda_values=None

    if type(lambda_input)==int: 
        if int(lambda_input)==3:
            output_lambda_values=[0.1127,0.5,0.88729]
        elif int(lambda_input)==5:
            output_lambda_values=[0.04691,0.23076,0.5,0.76923,0.95308]
        elif int(lambda_input)==7:
            output_lambda_values=[0.02544,0.12923,0.29707,0.5,0.70292,0.87076,0.97455]
        elif int(lambda_input)==9:
            output_lambda_values=[0.01592,0.08198,0.19331,0.33787,0.5,0.66213,0.80669,0.91802,0.98408]
        elif int(lambda_input)==12:
            output_lambda_values=[0.00922,0.04794,0.11505,0.20634,0.31608,0.43738,0.56262,0.68392,0.79366,0.88495,0.95206,0.99078]
        else:
            output_lambda_values=list(np.linspace(0,1,num=int(lambda_input)))
            output_lambda_values=[np.around(x,3) for x in output_lambda_values]

    elif type(lambda_input)==list:
        output_lambda_values=[np.around(float(x),3) for x in lambda_list]
        if output_lambda_values[-1]>1.0:
            raise Exception('Error: lambda windows cannot be above 1\n')

    output_lambda_values.sort(reverse=False)

    return output_lambda_values

def write_lambda_windows(media,protocol,schedule,ti_masks,ti_mask_len,hmass,equil_ns,prod_ns,monte_water=0,dir_1_name='core',dir_2_name='sec_lig'):

    prmtop_list=glob.glob('../%s*prmtop' % (media))

    # need the number of atoms in mol1 and mol2
    if protocol=='three-step':
        core=Chem.SDMolSupplier('../../'+dir_1_name+'/LIG.sdf',removeHs=False)[0]
        sec_lig=Chem.SDMolSupplier('../../'+dir_2_name+'/MOD.sdf',removeHs=False)[0]

    for prmtop in prmtop_list:
        prmtop_path=os.path.abspath(prmtop)

        sch=prmtop.split('/')[-1].split('.')[0]
        lambda_values=set_lambda_values(schedule[sch])

        os.mkdir(sch)
        os.chdir(sch)

        for i in range(0,len(lambda_values)):
            #print(i+1,lambda_values[i])

            write_ti_cluster_script(prmtop_path)

            with open('dir_list.dat','a') as f_out:
                f_out.write('%s\n' % ('lambda_'+str(i)))

            os.mkdir('lambda_'+str(i))
            os.chdir('lambda_'+str(i))

            if protocol=='three-step':
                write_ti_inputs(protocol,sch,ti_masks,ti_mask_len,lambda_values[i],lambda_values,hmass,equil_ns,prod_ns,monte_water,len(core.GetAtoms()),len(sec_lig.GetAtoms()))
            else:
                write_ti_inputs(protocol,sch,ti_masks,ti_mask_len,lambda_values[i],lambda_values,hmass,equil_ns,prod_ns,monte_water)

            # lambda window
            os.chdir('../')

        # schedule
        os.chdir('../')

def run_prod(df,protocol='one-step',ti_repeats=1,schedule={'complex_ligands':9,'solvent_ligands':9},hmass=True,equil_ns=1,prod_ns=5,monte_water=0):

    # submit all dir lig1->lig2
    for index,row in df.iterrows():
        name1=row['Name1']
        # rbfe
        if 'Name2' in list(df.columns):
            name2=row['Name2']
            pair_dir=name1+'~'+name2
        # absolute
        else:
            pair_dir=name1

        os.chdir(pair_dir)

        # Save TI masks
        if protocol not in ['absolute','absolute-three-step']:
            ti_masks=[]
            ti_mask_len=[]
            counter=1
            with open('TI_MASKS.dat','r') as f:
                for line in f:
                    ti_masks.append(':'+str(counter)+'@'+line.strip('\n'))
                    ti_mask_len.append(int(len(line.split(',')))-1)
                    counter+=1
        else:
            ti_masks=['','']
            ti_mask_len=[0,0]

        for media in ['complex','solvent']:
            os.chdir(media)

            for rep in range(1,ti_repeats+1):
                print('Submit production: %s %s repeat %d\n' % (pair_dir,media,rep))

                os.mkdir('%s_rep%d' % (protocol,rep))
                os.chdir('%s_rep%d' % (protocol,rep))

                # write lambda dir and submit
                write_lambda_windows(media,protocol,schedule,ti_masks,ti_mask_len,hmass,equil_ns,prod_ns,monte_water)

                # repeat
                os.chdir('../')

            # media
            os.chdir('../')
    
        # pair_dir
        os.chdir('../')
        
