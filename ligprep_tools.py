#!/usr/bin/env python
import sys
import os
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

    def query_atom_name(self,at_name):

        atom_out=None
        for atom in self._atoms:
            if atom.name==at_name:
                atom_out=atom
                break

        return atom_out

    def query_atom_type(self,type1):

        atom_out=None
        for atom in self._atoms:
            if atom.atom_type==type1:
                atom_out=atom
                break

        return atom_out

    def query_atom_index(self,atom):

        return self._atoms.index(atom)

    def query_atom_atomIdx(self,at_idx):
        try:
            self._atoms[at_idx]
            return at_idx
        except:
            return None

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

    def query_bond_atomIdx(self,at1_idx,at2_idx):

        bond_out=None
        for bond in self._bonds:
            if bond.atom1.idx==at1_idx and bond.atom2.idx==at2_idx:
                bond_out=bond
                break
            elif bond.atom1.idx==at2_idx and bond.atom2.idx==at1_idx:
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

    def query_angle_atomIdx(self,at1_idx,at2_idx,at3_idx):

        angle_out=None
        for angle in self._angles:
            if angle.atom1.idx==at1_idx and angle.atom2.idx==at2_idx and angle.atom3.idx==at3_idx:
                angle_out=angle
                break
            elif angle.atom1.idx==at3_idx and angle.atom2.idx==at2_idx and angle.atom3.idx==at1_idx:
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

    def query_dihedral_atomIdx(self,at1_idx,at2_idx,at3_idx,at4_idx):

        dihed_out=None
        for dihed in self._dihedrals:
            if dihed.atom1.idx==at1_idx and dihed.atom2.idx==at2_idx and dihed.atom3.idx==at3_idx and dihed.atom4.idx==at4_idx:
                dihed_out=dihed
                break
            elif dihed.atom1.idx==at4_idx and dihed.atom2.idx==at3_idx and dihed.atom3.idx==at2_idx and dihed.atom4.idx==at1_idx:
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
    def __init__(self,idx,name=None,atomic_num=None,atomic_weight=None,element=None,atom_type=None,atom_charge=None,lj_R=None,lj_E=None,hybrid=None,bond_count=None,x=None,y=None,z=None,ti_core=False):
        self.idx=idx
        self.atomic_num=atomic_num
        self.atomic_weight=atomic_weight
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
    def __init__(self,atom1,atom2,atom3,atom4,frc=None,period=None,phase=None,n_terms=None,improper=False):
        self.atom1=atom1
        self.atom2=atom2
        self.atom3=atom3
        self.atom4=atom4
        self.n_terms=n_terms
        self.improper=improper

        if frc is not None:
            self._frc=list(frc)
        else:
            self._frc=[]

        if period is not None:
            self._period=list(period)
        else:
            self._period=[]

        if phase is not None:
            self._phase=list(phase)
        else:
            self._phase=[]

    @property
    def frc(self):
        return self._frc

    @property
    def period(self):
        return self._period

    @property
    def phase(self):
        return self._phase

    def set_negative_period(self):
        for i in range(0,len(self._period)):
            if self._period[i]>1:
                self._period[i]=int(-1*abs(self._period[i]))

    def sort_dihed(self):
        self._period, self._frc, self._phase = zip(*sorted(zip(self._period, self._frc, self._phase)))
        self._period=list(self._period)
        self._frc=list(self._frc)
        self._phase=list(self._phase)

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

                if fields is not None:
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

# List all torsions present in a mol, skipping fragment
def enumerateSkipTorsions(mol,frag=None):
    torsionSmarts = '[!$(*#*)]~[!$(*#*)]~[!$(*#*)]~[!$(*#*)]'
    torsionQuery = Chem.MolFromSmarts(torsionSmarts)
    matches = mol.GetSubstructMatches(torsionQuery,uniquify=False)
    torsionList = []

    my_frag_tor=[]
    frag_skip=[]
    if frag is not None:
        if mol.HasSubstructMatch(frag):
            frag_idx=mol.GetSubstructMatch(frag)
            for val in frag_idx:
                frag_skip.append(val)

    for match in matches:
        if match not in torsionList:
            if not remove_frag_tor(mol,match,frag_skip):
                torsionList.append(match)
            else:
                my_frag_tor.append(match)

    return torsionList,my_frag_tor

# Skip a torsion if it is in the fragment we want to skip
def remove_frag_tor(mol,torsion,frag_skip):
    output=False

    tor_idx=[]
    for val in torsion:
        tor_idx.append(val)

    if (tor_idx[0] in frag_skip) and (tor_idx[1] in frag_skip) and (tor_idx[2] in frag_skip) and (tor_idx[3] in frag_skip):
        output=True
    # to avoid finding a hydrogen torsion
    elif (tor_idx[0] in frag_skip) and (tor_idx[1] in frag_skip) and (tor_idx[2] in frag_skip) and (int(mol.GetAtomWithIdx(tor_idx[3]).GetAtomicNum())==1):
        output=True
    elif (tor_idx[3] in frag_skip) and (tor_idx[2] in frag_skip) and (tor_idx[1] in frag_skip) and (int(mol.GetAtomWithIdx(tor_idx[0]).GetAtomicNum())==1):
        output=True
    return output

# unique torions
def unique_tor(torsionList,mol,mol_ff):
    exist_tor=[]
    exist_id=[]
    for torsion in torsionList:
        if check_torsion(mol,torsion):
            exist_add=str('%s-%s-%s-%s' % (mol_ff.atoms[torsion[0]].atom_type,mol_ff.atoms[torsion[1]].atom_type,mol_ff.atoms[torsion[2]].atom_type,mol_ff.atoms[torsion[3]].atom_type))
            reverse_add=str('%s-%s-%s-%s' % (mol_ff.atoms[torsion[3]].atom_type,mol_ff.atoms[torsion[2]].atom_type,mol_ff.atoms[torsion[1]].atom_type,mol_ff.atoms[torsion[0]].atom_type))
            if (exist_add not in exist_tor) and (reverse_add not in exist_tor):
                exist_tor.append(exist_add)
                exist_id.append(torsion)

    return exist_tor,exist_id

# Check if a torsion is not within a ring; is not of type H-C-X-X
def check_torsion(mol,torsion):
    central_type=mol.GetBondBetweenAtoms(torsion[1],torsion[2]).GetBondType()
    central_ring=(mol.GetAtomWithIdx(torsion[1]).IsInRing() and mol.GetAtomWithIdx(torsion[2]).IsInRing())

    atom_1=int(mol.GetAtomWithIdx(torsion[0]).GetAtomicNum())
    atom_2=int(mol.GetAtomWithIdx(torsion[1]).GetAtomicNum())
    atom_3=int(mol.GetAtomWithIdx(torsion[2]).GetAtomicNum())
    atom_4=int(mol.GetAtomWithIdx(torsion[3]).GetAtomicNum())

    a1_hybrid=str(mol.GetAtomWithIdx(torsion[1]).GetHybridization())
    b1_hybrid=str(mol.GetAtomWithIdx(torsion[2]).GetHybridization())

    result=False

    if (str(central_type)=='SINGLE' and central_ring==False):
        if (atom_1==1 and atom_2==6):
            result=False
        elif (atom_4==1 and atom_3==6):
            result=False
        else:
            result=True
    elif (str(central_type)=='SINGLE' and central_ring==True and idx_same_ring(mol,torsion)==False and idx_central_ring(mol,[torsion[1],torsion[2]])==False):
        if (a1_hybrid=='SP2') and (b1_hybrid=='SP2'):
            if (atom_1==1 and atom_2==6):
                result=False
            elif (atom_4==1 and atom_3==6):
                result=False
            else:
                result=True
        elif (a1_hybrid=='SP2') and (b1_hybrid=='SP3'):
            if (atom_1==1 and atom_2==6):
                result=False
            elif (atom_4==1 and atom_3==6):
                result=False
            else:
                result=True
        elif (a1_hybrid=='SP3') and (b1_hybrid=='SP2'):
            if (atom_1==1 and atom_2==6):
                result=False
            elif (atom_4==1 and atom_3==6):
                result=False
            else:
                result=True

        return result

def fix_duplicate(tor_list,mol):
    central_sort=[]
    central_unique=[]
    tor_mass=[]

    # central idx and tor mass
    for tor in tor_list:
        central_sort.append(sorted([tor.split('-')[1],tor.split('-')[2]],key=int))
        a=int(tor.split('-')[0])
        b=int(tor.split('-')[1])
        c=int(tor.split('-')[2])
        d=int(tor.split('-')[3])

        a_num=mol.GetAtomWithIdx(a-1).GetMass()
        b_num=mol.GetAtomWithIdx(b-1).GetMass()
        c_num=mol.GetAtomWithIdx(c-1).GetMass()
        d_num=mol.GetAtomWithIdx(d-1).GetMass()

        tor_mass.append(float(a_num+b_num+c_num+d_num))

    # unique central idx
    for tor in central_sort:
        if tor not in central_unique:
            central_unique.append(tor)

    per_tor=len(central_unique)*[None]

    counter=0
    for tor in central_unique:
        my_list=[]

        for i in range(0,len(tor_list)):
            if tor==central_sort[i]:
                my_list.append([tor_list[i],tor_mass[i]])

        per_tor[counter]=my_list
        counter+=1

    return per_tor

# check if 3 or more atoms are in the same ring
def idx_same_ring(mol,atom_idx):
    ring_data=[]
    for ring in mol.GetRingInfo().AtomRings():
        ring_data.append(ring)

    same_ring=False
    for i in range(0,len(mol.GetRingInfo().AtomRings())):
        found_in_ring=0
        for atom in ring_data[i]:
            for idx in atom_idx:
                if atom==idx:
                    found_in_ring+=1

            if found_in_ring>2:
                same_ring=True
                break

    return same_ring

def idx_central_ring(mol,atom_idx):
    ring_data=[]
    for ring in mol.GetRingInfo().AtomRings():
        ring_data.append(ring)

    same_ring=False
    for i in range(0,len(mol.GetRingInfo().AtomRings())):
        found_in_ring=0
        for atom in ring_data[i]:
            for idx in atom_idx:
                if atom==idx:
                    found_in_ring+=1

            if found_in_ring>1:
                same_ring=True
                break

    return same_ring


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
def build_parm(residue_name,ff,file_name,prmtop_name=None,frcmod_file=None):

    with open(file_name,'w') as f:
        if ff in ['gaff','gaff2','num']:
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
def save_atom(mol_ff,Atom_file,fields=None):

    with open(Atom_file,'r') as f:
        for line in f:
            if len(line.split())>10 and can_Int(line.split()[0]):
                at_idx=int(line.split()[0])-1

                name=line.split()[3]
                atom_type=line.split()[4]
                atom_charge=float(line.split()[9])
                atomic_weight=float(line.split()[8])
                atomic_num=int(line.split()[5])
                vdw_r=float(line.split()[6])
                vdw_e=float(line.split()[7])

                if mol_ff.query_atom_atomIdx(at_idx) is  None:
                    mol_ff.add_atom(Atom_ff(idx=at_idx))

                mol_ff.atoms[at_idx].lj_R=vdw_r
                mol_ff.atoms[at_idx].lj_E=vdw_e

                if mol_ff.atoms[at_idx].atomic_weight is None:
                    mol_ff.atoms[at_idx].atomic_weight=atomic_weight

                if mol_ff.atoms[at_idx].atomic_num is None:
                    mol_ff.atoms[at_idx].atomic_num=atomic_num

                if fields is not None:
                    if 'name' in fields:
                        mol_ff.atoms[at_idx].name=name
                    if 'type' in fields:
                        mol_ff.atoms[at_idx].atom_type=atom_type
                    if 'charge' in fields:
                        mol_ff.atoms[at_idx].atom_charge=atom_charge

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
                    mol_ff.add_bond(Bond_ff(mol_ff.atoms[at1_idx],mol_ff.atoms[at2_idx],frc_val,length_val))

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
                phase_val=round(float(line.split()[-3]),1)
    
                if line.split()[0]=='I':
                    improp=True

                if mol_ff.query_dihedral_atomIdx(at1_idx,at2_idx,at3_idx,at4_idx) is None:
                    mol_ff.add_dihedral(Dihedral_ff(mol_ff.atoms[at1_idx],mol_ff.atoms[at2_idx],mol_ff.atoms[at3_idx],mol_ff.atoms[at4_idx],frc=[frc_val],period=[period_val],phase=[phase_val],n_terms=1,improper=improp))
                else:
                    dihed_idx=mol_ff.query_dihedral_idx(mol_ff.query_dihedral_atomIdx(at1_idx,at2_idx,at3_idx,at4_idx))
                    if period_val not in mol_ff.dihedrals[dihed_idx].period:
                        if mol_ff.dihedrals[dihed_idx].n_terms==None:
                            mol_ff.dihedrals[dihed_idx].n_terms=1
                        else:
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

def norm_len(atom_char):
        if len(str(atom_char))>1:
                return str(atom_char)
        elif len(str(atom_char))==1:
                return str(atom_char)+' '

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
                        line=('%s-%s-%s-%s    1     %-6s %10s  %6s\n' % (norm_len(dihedrals[i].atom1.atom_type),norm_len(dihedrals[i].atom2.atom_type),norm_len(dihedrals[i].atom3.atom_type),norm_len(dihedrals[i].atom4.atom_type),float(dihedrals[i].frc[term]),float(dihedrals[i].phase[term]),float(dihedrals[i].period[term])))
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

def check_file(file_in):
    output=False
    if os.path.exists(file_in) and os.path.getsize(file_in)>0:
        output=True
    return output

def leaprc_out(leaprc_name,mol):
    hydrid={'SP3':'sp3','1':'sp3','SP3D':'sp3','SP2':'sp2','S':'sp3','SP':'sp'}
    periodic_lookup={'6':'C','1':'H','8':'O','7':'N','9':'F','15':'P','16':'S','17':'Cl','35':'Br','53':'I'}

    with open(leaprc_name,'w') as f:
        f.write('addAtomTypes {\n')
        for i in range(0,len(mol.GetAtoms())):
            f.write('{"%s" "%s" "%s"}\n' % (str(i+1),periodic_lookup[str(mol.GetAtomWithIdx(i).GetAtomicNum())],hydrid[str(mol.GetAtomWithIdx(i).GetHybridization())]))
        f.write('}\n')

##############################################################################

