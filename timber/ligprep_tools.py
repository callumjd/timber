# timber

import os
import numpy as np
import parmed as pmd
from rdkit import Chem
from rdkit.Chem import AllChem

##############################################################################

def check_file(file_in):
    output=False
    if os.path.exists(file_in) and os.path.getsize(file_in)>0:
        output=True
    return output

# write file to make amber .off
def make_off(mol_ff,file_name):

    with open(file_name,'w') as f:
        f.write('%s=loadpdb %s.pdb\n' % (mol_ff.name,mol_ff.name))

        for i in range(0,len(mol_ff.atoms)):
            f.write('set %s.1.%s type "%s"\n' % (mol_ff.name,mol_ff.atoms[i].name,mol_ff.atoms[i].atom_type))

        for i in range(0,len(mol_ff.atoms)):
            f.write('set %s.1.%s charge %s\n' % (mol_ff.name,mol_ff.atoms[i].name,mol_ff.atoms[i].atom_charge))

        for i in range(0,len(mol_ff.bonds)):
            f.write('bond %s.1.%s %s.1.%s\n' % (mol_ff.name,mol_ff.bonds[i].atom1.name,mol_ff.name,mol_ff.bonds[i].atom2.name))

        f.write('saveoff %s %s.off\n' % (mol_ff.name,mol_ff.name))
        f.write('quit\n')

# run antechamber
# if net_charge is defined, calculate am1-bcc charges
def run_antechamber(input_file,residue_name='UNL',ff='gaff2',net_charge=None):

    file_format=input_file.split('.')[-1]

    try:
        assert(file_format in ['mol2','sdf'])
    except:
        raise Exception('File format error')

    try:
        assert(ff in ['gaff','gaff2'])
    except:
        raise Exception('Force field undefined')

    if net_charge==None:
        os.system('antechamber -i %s -fi %s -o %s.mol2 -fo mol2 -rn %s -at %s -s 0 -pf y -dr no' % (input_file,file_format,residue_name,residue_name,ff))
    else:
        os.system('antechamber -i %s -fi %s -o %s.mol2 -fo mol2 -rn %s -nc %d -c bcc -at %s -s 0 -pf y -dr no' % (input_file,file_format,residue_name,residue_name,net_charge,ff))
        os.system('rm sqm.in sqm.out sqm.pdb')

    os.system('parmchk2 -i %s.mol2 -f mol2 -o missing_%s.frcmod -s %s' % (residue_name,ff,ff))

    # leap will fail if file does not exist
    if not check_file('missing_%s.frcmod' % (ff)):
        os.system('missing_%s.frcmod' % (ff))

    # clean SDF file for rdkit
    os.system('antechamber -i %s.mol2 -fi mol2 -o %s.sdf -fo sdf -s 0 -pf y -dr no' % (residue_name,residue_name))

# write file to build amber prmtop
def build_parm(residue_name='UNL',ff='gaff2',file_name='make_lig.leap',prmtop_name=None,frcmod_file=None):

    with open(file_name,'w') as f:
        if ff in ['gaff','gaff2']:
            f.write('source leaprc.%s\n' % (ff))
        else:
            f.write('source leaprc.protein.ff14SB\n')
        if frcmod_file is not None:
            if type(frcmod_file)==str:
                frcmod_file=[frcmod_file]
            for frc_file in frcmod_file:
                f.write('loadamberparams %s\n' % (frc_file))
        f.write('loadoff %s.off\n' % (residue_name))
        f.write('mol=loadpdb %s.pdb\n' % (residue_name))
        if prmtop_name is not None:
            f.write('saveamberparm mol %s.prmtop inpcrd\n' % (prmtop_name))
        else:
            f.write('saveamberparm mol prmtop inpcrd\n')
        f.write('quit\n')

# Take .mol2 information to mol_ff
def Info_Mol2(mol2_file,mol_ff,n_atoms,fields=None):

    start=0
    with open(mol2_file,'r') as f:
        for line in f:
            start+=1
            if '<TRIPOS>ATOM' in line:
                break

    counter=0
    with open(mol2_file,'r') as f:
        for line in f:
            if (start-1<counter<n_atoms+start):
                name=line.split()[1]
                atom_type=line.split()[5]
                atom_charge=line.split()[8]

                if fields is not None:
                    if 'name' in fields:
                        mol_ff.atoms[counter-start].name=str(name)
                    if 'type' in fields:
                        mol_ff.atoms[counter-start].atom_type=str(atom_type)
                    if 'charge' in fields:
                        mol_ff.atoms[counter-start].atom_charge=float(atom_charge)
            counter+=1

    return mol_ff

# Take .off information to mol_ff
def Info_OFF(off_file,mol_ff,n_atoms,fields=None):

    start=0
    with open(off_file,'r') as f:
        for line in f:
            start+=1
            if '!entry' in line:
                break

    counter=0
    with open(off_file,'r') as f:
        for line in f:
            if (start-1<counter<n_atoms+start):
                name=line.split()[0].strip('"')
                atom_type=line.split()[1].strip('"')
                atom_charge=float(line.split()[7])

                if fields is not None:
                    if 'name' in fields:
                        mol_ff.atoms[counter-start].name=str(name)
                    if 'type' in fields:
                        mol_ff.atoms[counter-start].atom_type=str(atom_type)
                    if 'charge' in fields:
                        mol_ff.atoms[counter-start].atom_charge=float(atom_charge)
            counter+=1

    return mol_ff

def write_rd_pdb(mol_ff,rd_mol,residue_name,output_file):

    counter=0
    for atom in rd_mol.GetAtoms():
        mi = Chem.AtomPDBResidueInfo()
        mi.SetName(mol_ff.atoms[counter].name)
        # the rdkit PDB residue name has incorrect whitespace
        mi.SetResidueName(''.ljust(4-len(mol_ff.atoms[counter].name))+residue_name)
        mi.SetResidueNumber(1)
        mi.SetIsHeteroAtom(False)
        atom.SetMonomerInfo(mi)

        counter+=1

    Chem.MolToPDBFile(rd_mol,output_file,flavor=2)

    # CONECT records break leap
    pdb_data=[]
    with open(output_file,'r') as f:
        for line in f:
            if line.split()[0]=='ATOM' or line.split()[0]=='HETATM':
                pdb_data.append(line)

    with open(output_file,'w') as f:
        for line in pdb_data:
            f.write(line)

def write_frcmod(file_name,types,bonds,angles,dihedrals):
    with open(file_name,'w') as f:
        f.write('missing\n')
        f.write('MASS\n')
        if types is not None:
            for i in range(0,len(types)):
                line=('%s    %6s\n' % (norm_len(types[i].atom_type),float(types[i].atomic_weight)))
                f.write(line)

        f.write('\n')
        f.write('BOND\n')
        if bonds is not None:
            for i in range(0,len(bonds)):
                line=('%s-%s    %s    %s\n' % (norm_len(bonds[i].atom1.atom_type),norm_len(bonds[i].atom2.atom_type),bonds[i].frc,bonds[i].length))
                f.write(line)

        f.write('\n')
        f.write('ANGLE\n')
        if angles is not None:
            for i in range(0,len(angles)):
                line=('%s-%s-%s    %s    %s\n' % (norm_len(angles[i].atom1.atom_type),norm_len(angles[i].atom2.atom_type),norm_len(angles[i].atom3.atom_type),angles[i].frc,angles[i].ref_angle))
                f.write(line)

        f.write('\n')
        f.write('DIHE\n')
        if dihedrals is not None:
            for i in range(0,len(dihedrals)):
                if dihedrals[i].improper==False:
                    for term in range(0,dihedrals[i].n_terms):
                        line=('%s-%s-%s-%s    %d     %-6s %10s  %6s\n' % (norm_len(dihedrals[i].atom1.atom_type),norm_len(dihedrals[i].atom2.atom_type),norm_len(dihedrals[i].atom3.atom_type),norm_len(dihedrals[i].atom4.atom_type),int(dihedrals[i].idiv),float(dihedrals[i].frc[term]),float(dihedrals[i].phase[term]),float(dihedrals[i].period[term])))
                        f.write(line)

        f.write('\n')
        f.write('IMPROPER\n')
        if dihedrals is not None:
            for i in range(0,len(dihedrals)):
                if dihedrals[i].improper==True:
                    for term in range(0,dihedrals[i].n_terms):
                        line=('%s-%s-%s-%s    %6s %10s  %6s\n' % (norm_len(dihedrals[i].atom1.atom_type),norm_len(dihedrals[i].atom2.atom_type),norm_len(dihedrals[i].atom3.atom_type),norm_len(dihedrals[i].atom4.atom_type),float(dihedrals[i].frc[term]),float(dihedrals[i].phase[term]),float(dihedrals[i].period[term])))
                        f.write(line)

        f.write('\n')
        f.write('NONB\n')
        if types is not None:
            for i in range(0,len(types)):
                line=('%s    %-6s    %-6s\n' % (norm_len(types[i].atom_type),float(types[i].lj_R),float(types[i].lj_E)))
                f.write(line)
        f.write('\n')

def leaprc_out(leaprc_name,mol,mol_ff):
    hydrid={'SP3':'sp3','1':'sp3','SP3D':'sp3','SP2':'sp2','S':'sp3','SP':'sp'}
    periodic_lookup={'6':'C','1':'H','8':'O','7':'N','9':'F','15':'P','16':'S','17':'Cl','35':'Br','53':'I'}

    with open(leaprc_name,'w') as f:
        f.write('addAtomTypes {\n')
        for i in range(0,len(mol.GetAtoms())):
            f.write('{"%s" "%s" "%s"}\n' % (mol_ff.atoms[i].atom_type,periodic_lookup[str(mol.GetAtomWithIdx(i).GetAtomicNum())],hydrid[str(mol.GetAtomWithIdx(i).GetHybridization())]))
        f.write('}\n')

# update rdmol coords from pdb file
def update_mol_coords_pdb(rdmol,pdb_file):

    copy_mol=copy.deepcopy(rdmol)
    # get updated MD coords for ligand
    lig_coords=[]
    with open(pdb_file,'r') as f:
        for line in f:
            if (len(line.split())>1 and (line.split()[0]=='ATOM')):
                idx=int(line.split()[1])
                x=float(line.split()[5])
                y=float(line.split()[6])
                z=float(line.split()[7])

                copy_mol.GetConformer().SetAtomPosition(idx-1,(x,y,z))
                lig_coords.append(Coord(x,y,z))

    return copy_mol,lig_coords

# run hydrogen mass repartitioning
def setup_hmass(prmtop):

    os.system('mv %s nohmass.%s' % (prmtop,prmtop))    

    my_prmtop=pmd.amber.AmberParm('nohmass.'+prmtop)

    action=pmd.tools.HMassRepartition(my_prmtop)
    action.execute()
    my_prmtop.write_parm(prmtop)

