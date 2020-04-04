#!/usr/bin/env python
import sys
import os
import copy
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

##############################################################################

class Molecule_ff():

    def __init__(self,name=None):
        if name is None:
            name = ''
        self._name=name

        self._atoms=list()
        self._bonds=list()
        self._angles=list()
        self._dihedrals=list()

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self,value):
        self._name=value

    @property
    def atoms(self):
        return self._atoms

    @property
    def bonds(self):
        return self._bonds

    @property
    def angles(self):
        return self._angles

    @property
    def dihedrals(self):
        return self._dihedrals

    def add_atom(self,atom):
        self._atoms.append(atom)

    def add_bond(self,bond):
        if bond not in self._bonds:
            self._bonds.append(bond)

    def add_angle(self,angle):
        if angle not in self._angles:
            self._angles.append(angle)

    def add_dihedral(self,dihedral):
        if dihedral not in self._dihedrals:
            self._dihedrals.append(dihedral)

    def query_atom_type(self,type1):

        atom_out=None
        for atom in self._atoms:
            if atom.atom_type==type1:
                atom_out=atom
                break

        return atom_out

    def query_atom_index(self,atom):

        return self._atoms.index(atom)

    def query_bond_type(self,type1,type2):

        bond_out=None
        for bond in self._bonds:
            if bond.atom1.atom_type==type1 and bond.atom2.atom_type==type2:
                bond_out=bond
                break
            elif bond.atom1.atom_type==type2 and bond.atom2.atom_type==type1:
                bond_out=bond
                break

        return bond_out

    def query_bond_atomIdx(self,at_idx1,at_idx2):

        bond_out=None
        for bond in self._bonds:
            if bond.atom1.idx==at_idx1 and bond.atom2.idx==at_idx2:
                bond_out=bond
                break
            elif bond.atom1.idx==at_idx2 and bond.atom2.idx==at_idx1:
                bond_out=bond
                break

        return bond_out

    def query_bond_index(self,bond):

        return self._bonds.index(bond)

    def query_angle_type(self,type1,type2,type3):

        angle_out=None
        for angle in self._angles:
            if angle.atom1.atom_type==type1 and angle.atom2.atom_type==type2 and angle.atom3.atom_type==type3:
                angle_out=angle
                break
            elif angle.atom1.atom_type==type3 and angle.atom2.atom_type==type2 and angle.atom3.atom_type==type1:
                angle_out=angle
                break

        return angle_out

    def query_angle_atomIdx(self,at_idx1,at_idx2,at_idx3):

        angle_out=None
        for angle in self._angles:
            if angle.atom1.idx==at_idx1 and angle.atom2.idx==at_idx2 and angle.atom3.idx==at_idx3:
                angle_out=angle
                break
            elif angle.atom1.idx==at_idx3 and angle.atom2.idx==at_idx2 and angle.atom3.idx==at_idx1:
                angle_out=angle
                break

        return angle_out

    def query_angle_index(self,angle):

        return self._angles.index(angle)

    def query_dihedral_type(self,type1,type2,type3,type4):

        dihed_out=None
        for dihed in self._dihedrals:
            if dihed.atom1.atom_type==type1 and dihed.atom2.atom_type==type2 and dihed.atom3.atom_type==type3 and dihed.atom4.atom_type==type4:
                dihed_out=dihed
                break
            elif dihed.atom1.atom_type==type4 and dihed.atom2.atom_type==type3 and dihed.atom3.atom_type==type2 and dihed.atom4.atom_type==type1:
                dihed_out=dihed
                break

        return dihed_out

    def query_dihedral_atomIdx(self,at_idx1,at_idx2,at_idx3,at_idx4):

        dihed_out=None
        for dihed in self._dihedrals:
            if dihed.atom1.idx==at_idx1 and dihed.atom2.idx==at_idx2 and dihed.atom3.idx==at_idx3 and dihed.atom4.idx==at_idx4:
                dihed_out=dihed
                break
            elif dihed.atom1.idx==at_idx4 and dihed.atom2.idx==at_idx3 and dihed.atom3.idx==at_idx2 and dihed.atom4.idx==at_idx1:
                dihed_out=dihed
                break

        return dihed_out

    def query_dihedral_idx(self,dihed):

        return self._dihedrals.index(dihed)

    def clear_atoms(self):
        self._atoms=list()

    def clear_bonds(self):
        self._bonds=list()

    def clear_angles(self):
        self._angles=list()

    def clear_dihedrals(self):
        self._dihedrals=list()

    def reorder_atoms(self,new_order):
        self._atoms=[self._atoms[i] for i in new_order]

class Atom_ff(object):
    def __init__(self,idx,name=None,mass=None,element=None,atom_type=None,atom_charge=None,lj_R=None,lj_E=None,hybrid=None,bond_count=None,x=None,y=None,z=None,ti_core=None):
        self.idx=idx
        self.name=name
        self.mass=mass
        self.element=element
        self.atom_type=atom_type
        self.atom_charge=atom_charge
        self.lj_R=lj_R
        self.lj_E=lj_E
        self.hybrid=hybrid
        self.bond_count=bond_count
        self.x=x
        self.y=y
        self.z=z
        self.ti_core=ti_core

class Bond_ff(object):
    def __init__(self,atom1,atom2,frc=None,length=None):
        self.atom1=atom1
        self.atom2=atom2
        self.frc=frc
        self.length=length

class Angle_ff(object):
    def __init__(self,atom1,atom2,atom3,frc=None,ref_angle=None):
        self.atom1=atom1
        self.atom2=atom2
        self.atom3=atom3
        self.frc=frc
        self.ref_angle=ref_angle

class Dihedral_ff(object):
    def __init__(self,atom1,atom2,atom3,atom4,frc=None,period=None,phase=None,n_terms=None,improper=None):
        self.atom1=atom1
        self.atom2=atom2
        self.atom3=atom3
        self.atom4=atom4
        self.n_terms=n_terms
        self.improper=improper

        if frc is not None:
            self._frc=list(frc)
        if period is not None:
            self._period=list(period)
        if phase is not None:
            self._phase=list(phase)

    @property
    def frc(self):
        return self._frc

    @property
    def period(self):
        return self._period

    @property
    def phase(self):
        return self._phase

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

                if 'name' in fields:
                    mol_ff.atoms[counter-start].name=name
                if 'type' in fields:
                    mol_ff.atoms[counter-start].atom_type=atom_type
                if 'charge' in fields:
                    mol_ff.atoms[counter-start].atom_charge=atom_charge
            counter+=1

    return mol_ff

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

                if 'name' in fields:
                    mol_ff.atoms[counter-start].name=name
                if 'type' in fields:
                    mol_ff.atoms[counter-start].atom_type=atom_type
                if 'charge' in fields:
                    mol_ff.atoms[counter-start].atom_charge=atom_charge
            counter+=1

    return mol_ff

# add bonds to Molecule_FF
def add_rd_bonds(rd_mol,mol_ff):

    unique_bonds=[]
    for rd_atom in rd_mol.GetAtoms():
        for rd_bond in rd_atom.GetBonds():
            a1=rd_bond.GetBeginAtomIdx()
            b1=rd_bond.GetEndAtomIdx()

            if a1<b1:
                if (a1,b1) not in unique_bonds:
                    unique_bonds.append((a1,b1))
                    mol_ff.add_bond(Bond_ff(mol_ff.atoms[a1],mol_ff.atoms[b1]))
            elif b1<a1:
                if (b1,a1) not in unique_bonds:
                    unique_bonds.append((b1,a1))
                    mol_ff.add_bond(Bond_ff(mol_ff.atoms[b1],mol_ff.atoms[a1]))

    return mol_ff

# add angles to Molecule_FF
def add_rd_angles(rd_mol,mol_ff):

    unique_angles=[]
    angleList=enumerateAngles(rd_mol)

    for ang in angleList:
        if ang[0]<ang[2]:
            if (ang[0],ang[1],ang[2]) not in unique_angles:
                unique_angles.append((ang[0],ang[1],ang[2]))
                mol_ff.add_angle(Angle_ff(mol_ff.atoms[ang[0]],mol_ff.atoms[ang[1]],mol_ff.atoms[ang[2]]))
            else:
                if (ang[2],ang[1],ang[0]) not in unique_angles:
                    unique_angles.append((ang[2],ang[1],ang[0]))
                    mol_ff.add_angle(Angle_ff(mol_ff.atoms[ang[2]],mol_ff.atoms[ang[1]],mol_ff.atoms[ang[0]]))

    return mol_ff

# add dihedrals to Molecule_FF
def add_rd_torsions(rd_mol,mol_ff):

    unique_torsions=[]
    torsionList=enumerateTorsions(rd_mol)

    for tor in torsionList:
        if tor[0]<tor[3]:
            if (tor[0],tor[1],tor[2],tor[3]) not in unique_torsions:
                unique_torsions.append((tor[0],tor[1],tor[2],tor[3]))
                mol_ff.add_dihedral(Dihedral_ff(mol_ff.atoms[tor[0]],mol_ff.atoms[tor[1]],mol_ff.atoms[tor[2]],mol_ff.atoms[tor[3]]))
        else:
            if (tor[3],tor[2],tor[1],tor[0]) not in unique_torsions:
                unique_torsions.append((tor[3],tor[2],tor[1],tor[0]))
                mol_ff.add_dihedral(Dihedral_ff(mol_ff.atoms[tor[3]],mol_ff.atoms[tor[2]],mol_ff.atoms[tor[1]],mol_ff.atoms[tor[0]]))

    return mol_ff

# List all angles present in a rd_mol object
def enumerateAngles(mol):
        angleSmarts = '[*]~[*]~[*]'
        angleQuery = Chem.MolFromSmarts(angleSmarts)
        matches = mol.GetSubstructMatches(angleQuery)
        angleList = []

        for match in matches:
                angleList.append(match)

        return angleList

# List all torsions present in a rd_mol object
def enumerateTorsions(mol):
    torsionSmarts = '[*]~[*]~[*]~[*]'
    torsionQuery = Chem.MolFromSmarts(torsionSmarts)
    matches = mol.GetSubstructMatches(torsionQuery)
    torsionList = []

    for match in matches:
        torsionList.append(match)

    return torsionList

# write file to make amber .off
def make_off(mol_ff,file_name):

    with open(file_name,'w') as f:
        f.write('%s=loadpdb %s.pdb\n' % (mol_ff.name,mol_ff.name))

        for i in range(0,len(mol_ff.atoms)):
            f.write('set %s.1.%s type %s\n' % (mol_ff.name,mol_ff.atoms[i].name,mol_ff.atoms[i].atom_type))

        for i in range(0,len(mol_ff.atoms)):
            f.write('set %s.1.%s charge %s\n' % (mol_ff.name,mol_ff.atoms[i].name,mol_ff.atoms[i].atom_charge))

        for i in range(0,len(mol_ff.bonds)):
            f.write('bond %s.1.%s %s.1.%s\n' % (mol_ff.name,mol_ff.bonds[i].atom1.name,mol_ff.name,mol_ff.bonds[i].atom2.name))

        f.write('saveoff %s %s.off\n' % (mol_ff.name,mol_ff.name))
        f.write('quit\n')

# run antechamber
def run_antechamber(input_file,residue_name,ff,net_charge=None,clean_sdf=False):

    file_format=input_file.split('.')[1]
    if file_format in ['mol2','sdf'] and ff in ['gaff','gaff2']:
        if net_charge==None:
            os.system('antechamber -i %s -fi %s -o %s.mol2 -fo mol2 -rn %s -at %s -s 0 -pf y' % (input_file,file_format,residue_name,residue_name,ff))
        else:
            os.system('antechamber -i %s -fi %s -o %s.mol2 -fo mol2 -rn %s -nc %d -c bcc -at %s -s 0 -pf y' % (input_file,file_format,residue_name,residue_name,net_charge,ff))
        os.system('parmchk2 -i %s.mol2 -f mol2 -o missing_gaff.frcmod -at %s' % (residue_name,ff))

    # clean SDF file for rdkit
    if clean_sdf:
        os.system('antechamber -i %s.mol2 -fi mol2 -o %s.sdf -fo sdf -s 0 -pf y' % (residue_name,residue_name))

# write file to build amber prmtop
def build_parm(residue_name,ff,file_name):
    if ff=='gaff':
        with open(file_name,'w') as f:
            f.write('source leaprc.gaff\n')
            f.write('loadamberparams missing_gaff.frcmod\n')
            f.write('loadoff %s.off\n' % (residue_name))
            f.write('mol=loadpdb %s.pdb\n' % (residue_name))
            f.write('saveamberparm mol prmtop inpcrd\n')
            f.write('quit\n')

# parmed info
def get_info(file_name,which):

    with open(file_name,'w') as f:
        f.write('parm prmtop\n')
        f.write('print%s *\n' % (which))
        f.write('quit\n')

def can_Int(val):
    try:
        int(val)
        return True
    except:
        return False

# save atom info from parmed output
def save_atom(mol_ff,Atom_file):

    with open(Atom_file,'r') as f:
        for line in f:
            if len(line.split())>10 and can_Int(line.split()[0]):
                at_idx=int(line.split()[0])-1

                vdw_r=float(line.split()[6])
                vdw_e=float(line.split()[7])
                at_mass=float(line.split()[8])

                mol_ff.atoms[at_idx].lj_R=vdw_r
                mol_ff.atoms[at_idx].lj_E=vdw_e
                mol_ff.atoms[at_idx].mass=at_mass

    return mol_ff

# save bond info from parmed output
def save_bond(mol_ff,Bond_file):

    with open(Bond_file,'r') as f:
        for line in f:
            if len(line.split())>9 and can_Int(line.split()[0]):
                at1_idx=int(line.split()[0])-1
                at2_idx=int(line.split()[4])-1

                frc_val=float(line.split()[9])
                length_val=float(line.split()[8])

                if mol_ff.query_bond_atomIdx(at1_idx,at2_idx):
                    bond_idx=mol_ff.query_bond_index(mol_ff.query_bond_atomIdx(at1_idx,at2_idx))
                    mol_ff.bonds[bond_idx].frc=frc_val
                    mol_ff.bonds[bond_idx].length=length_val

                else:
                    mol_ff.add_bond(mol_ff.atoms[at1_idx],mol_ff.atoms[at2_idx],frc_val,length_val)

    return mol_ff

# save angle info from parmed output
def save_angle(mol_ff,Angle_file):

    with open(Angle_file,'r') as f:
        for line in f:
            if len(line.split())>12 and can_Int(line.split()[0]):
                at1_idx=int(line.split()[0])-1
                at2_idx=int(line.split()[4])-1
                at3_idx=int(line.split()[8])-1

                frc_val=float(line.split()[12])
                ref_angle_val=float(line.split()[13])

                if mol_ff.query_angle_atomIdx(at1_idx,at2_idx,at3_idx):
                    angle_idx=mol_ff.query_angle_index(mol_ff.query_angle_atomIdx(at1_idx,at2_idx,at3_idx))
                    mol_ff.angles[angle_idx].frc=frc_val
                    mol_ff.angles[angle_idx].ref_angle=ref_angle_val

                else:
                    mol_ff.add_angle(Angle_ff(mol_ff.atoms[at1_idx],mol_ff.atoms[at2_idx],mol_ff.atoms[at3_idx],frc_val,ref_angle_val))

    return mol_ff

# save dihedral info from parmed output
def save_dihed(mol_ff,Dihed_file):

    with open(Dihed_file,'r') as f:
        for line in f:
            if len(line.split())>20:
                improp=False
                at1_idx=int(line.split()[-21])-1
                at2_idx=int(line.split()[-17])-1
                at3_idx=int(line.split()[-13])-1
                at4_idx=int(line.split()[-9])-1

                frc_val=float(line.split()[-5])
                period_val=int(line.split()[-4].split('.')[0])
                phase_val=int(line.split()[-3].split('.')[0])
    
                if line.split()[0]=='I':
                    improp=True

                if mol_ff.query_dihedral_atomIdx(at1_idx,at2_idx,at3_idx,at4_idx) is None:
                    mol_ff.add_dihedral(Dihedral_ff(mol_ff.atoms[at1_idx],mol_ff.atoms[at2_idx],mol_ff.atoms[at3_idx],mol_ff.atoms[at4_idx],frc=[frc_val],period=[period_val],phase=[phase_val],n_terms=1,improper=improp))
                else:
                    dihed_idx=mol_ff.query_dihedral_idx(mol_ff.query_dihedral_atomIdx(at1_idx,at2_idx,at3_idx,at4_idx))
                    if period_val not in mol_ff.dihedrals[dihed_idx].period:
                        mol_ff.dihedrals[dihed_idx].n_terms=mol_ff.dihedrals[dihed_idx].n_terms+1
                        mol_ff.dihedrals[dihed_idx].frc.append(frc_val)
                        mol_ff.dihedrals[dihed_idx].period.append(period_val)
                        mol_ff.dihedrals[dihed_idx].phase.append(phase_val)

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

##############################################################################

