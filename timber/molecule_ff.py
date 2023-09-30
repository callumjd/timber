# timber

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem,rdMolTransforms,rdmolops
import parmed as pmd
from .ligprep_tools import check_file
from .geometry import Coord

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

    def query_13angle_atomIdx(self,at1_idx,at3_idx):

        angle_out=None
        for angle in self._angles:
            if angle.atom1.idx==at1_idx and angle.atom3.idx==at3_idx:
                angle_out=angle
                break
            elif angle.atom1.idx==at3_idx and angle.atom3.idx==at1_idx:
                angle_out=angle
                break

        return angle_out

    def query_any_angle_atomIdx(self,at1_idx):

        angle_out=[]
        for angle in self._angles:
            if angle.atom1.idx==at1_idx:
                angle_out.append(angle)
            elif angle.atom2.idx==at1_idx:
                angle_out.append(angle)
            elif angle.atom3.idx==at1_idx:
                angle_out.append(angle)

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

    def query_nb14_atomIdx(self,at1_idx,at4_idx):

        dihed_out=None
        for dihed in self._dihedrals:
            if dihed.atom1.idx==at1_idx and dihed.atom4.idx==at4_idx:
                dihed_out=dihed
                break
            elif dihed.atom1.idx==at4_idx and dihed.atom4.idx==at1_idx:
                dihed_out=dihed
                break

        return dihed_out

    def query_any_dihedral_atomIdx(self,at1_idx):

        dihed_out=[]
        for dihed in self._dihedrals:
            if dihed.atom1.idx==at1_idx:
                dihed_out.append(dihed)
            elif dihed.atom2.idx==at1_idx:
                dihed_out.append(dihed)
            elif dihed.atom3.idx==at1_idx:
                dihed_out.append(dihed)
            elif dihed.atom4.idx==at1_idx:
                dihed_out.append(dihed)

        return dihed_out

    def query_dihedral_idx(self,dihed):

        return self._dihedrals.index(dihed)

    # given a dihedral, assign multiplicity as required for frcmod
    def assign_dihedral_idiv(self,dihed):

        idx=0

        b=dihed.atom2.idx
        c=dihed.atom3.idx

        t1=dihed.atom1.atom_type
        t2=dihed.atom2.atom_type
        t3=dihed.atom3.atom_type
        t4=dihed.atom4.atom_type

        for mydihed in self._dihedrals:
            if (mydihed.atom2.idx==b) and (mydihed.atom3.idx==c):
                if (mydihed.atom1.atom_type==t1) and (mydihed.atom2.atom_type==t2) and (mydihed.atom3.atom_type==t3) and (mydihed.atom4.atom_type==t4):
                    idx+=1
            if (mydihed.atom3.idx==b) and (mydihed.atom2.idx==c):
                if (mydihed.atom4.atom_type==t1) and (mydihed.atom3.atom_type==t2) and (mydihed.atom2.atom_type==t3) and (mydihed.atom1.atom_type==t4):
                    idx+=1

        dihed.idiv=idx

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
    def __init__(self,idx,name=None,atomic_num=None,atomic_weight=None,element=None,atom_type=None,atom_charge=None,lj_R=None,lj_E=None,hybrid=None,bond_count=None,x=None,y=None,z=None,ti_core=False,ti_sc=False,resi=None):
        self.idx=idx
        self.name=name
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
        self.ti_sc=ti_sc
        self.resi=resi
   
    @property
    def xyz(self):
        return Coord(self.x,self.y,self.z)

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
    def __init__(self,atom1,atom2,atom3,atom4,idiv=1,frc=None,period=None,phase=None,n_terms=None,improper=False):
        self.atom1=atom1
        self.atom2=atom2
        self.atom3=atom3
        self.atom4=atom4
        self.idiv=idiv
        self.improper=improper

        if n_terms:
            self._n_terms=n_terms
        else:
            self._n_terms=0

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

    @property
    def n_terms(self):
        return self._n_terms

    @n_terms.setter
    def n_terms(self,value):
        self._n_terms=value

    @property
    def idiv(self):
        return self._idiv

    @idiv.setter
    def idiv(self,value):
        self._idiv=value

    def set_negative_period(self):
        for i in range(0,len(self._period)):
            if self._period[i]>min(self._period):
                self._period[i]=int(-1*abs(self._period[i]))

    def sort_dihed(self):
        self._period, self._frc, self._phase = zip(*sorted(zip(self._period, self._frc, self._phase)))
        self._period=list(self._period)
        self._frc=list(self._frc)
        self._phase=list(self._phase)

    def update_dihed_frc(self,frc_vals):
        assert len(frc_vals)==self._n_terms
        for i in range(0,len(frc_vals)):
            self._frc[i]=frc_vals[i]

    def update_dihed_period(self,period_vals):
        assert len(period_vals)==self._n_terms
        for i in range(0,len(period_vals)):
            self._period[i]=period_vals[i]

    def update_dihed_phase(self,phase_vals):
        assert len(phase_vals)==self._n_terms
        for i in range(0,len(phase_vals)):
            self._phase[i]=phase_vals[i]

    def zero_dihed(self,torIdx):
        self._period=[-5,-4,-3,-2,1]
        self._frc=[0,0,0,0,0]
        self._phase=[0,0,0,0,0]
        self._n_terms=5

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

def get_rdkit_info(rd_mol,mol_ff):

    for at in rd_mol.GetAtoms():
        pos=rd_mol.GetConformer().GetAtomPosition(at.GetIdx())
        mol_ff.add_atom(Atom_ff(idx=int(at.GetIdx()),atomic_num=int(at.GetAtomicNum()),atomic_weight=float(Chem.GetPeriodicTable().GetAtomicWeight(at.GetAtomicNum())),hybrid=at.GetHybridization(),bond_count=len(at.GetBonds()),x=float(pos.x),y=float(pos.y),z=float(pos.z)))
    mol_ff=add_rd_bonds(rd_mol,mol_ff)
    mol_ff=add_rd_angles(rd_mol,mol_ff)
    mol_ff=add_rd_torsions(rd_mol,mol_ff)

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

# save parmed data to mol_ff
def get_parmed_info(mol_ff,prmtop='prmtop',parmed_info=['atoms','bonds','angles','dihedrals']):

    if not check_file(prmtop):
        raise Exception('Missing file')

    mol_pmd=pmd.load_file(prmtop)

    if 'atoms' in parmed_info:
        mol_ff=save_parmed_atoms(mol_ff,mol_pmd)
    if 'bonds' in parmed_info:
        mol_ff=save_parmed_bonds(mol_ff,mol_pmd)
    if 'angles' in parmed_info:
        mol_ff=save_parmed_angles(mol_ff,mol_pmd)
    if 'dihedrals' in parmed_info:
        mol_ff=save_parmed_dihedrals(mol_ff,mol_pmd)

    return mol_ff

def save_parmed_atoms(mol_ff,mol_pmd):

    for at_idx in range(0,len(mol_pmd.atoms)):

        if mol_ff.query_atom_atomIdx(at_idx) is None:
            mol_ff.add_atom(Atom_ff(idx=at_idx))

        mol_ff.atoms[at_idx].name=str(mol_pmd.atoms[at_idx].name)
        mol_ff.atoms[at_idx].atom_type=str(mol_pmd.atoms[at_idx].type)
        mol_ff.atoms[at_idx].atom_charge=float(mol_pmd.atoms[at_idx].charge)
        mol_ff.atoms[at_idx].lj_R=float(mol_pmd.atoms[at_idx].rmin)
        mol_ff.atoms[at_idx].lj_E=float(mol_pmd.atoms[at_idx].epsilon)
        mol_ff.atoms[at_idx].resi=int(mol_pmd.atoms[at_idx].residue.number)+1
        mol_ff.atoms[at_idx].atomic_weight=float(mol_pmd.atoms[at_idx].mass)
        mol_ff.atoms[at_idx].atomic_num=int(mol_pmd.atoms[at_idx].atomic_number)

    return mol_ff

def save_parmed_bonds(mol_ff,mol_pmd):

    for bond_pmd in mol_pmd.bonds:
        at1_idx=bond_pmd.atom1.idx
        at2_idx=bond_pmd.atom2.idx
        frc_val=float(bond_pmd.type.k)
        length_val=float(bond_pmd.type.req)

        if mol_ff.query_bond_atomIdx(at1_idx,at2_idx):
            bond_idx=mol_ff.query_bond_index(mol_ff.query_bond_atomIdx(at1_idx,at2_idx))
            mol_ff.bonds[bond_idx].frc=frc_val
            mol_ff.bonds[bond_idx].length=length_val

        else:
            mol_ff.add_bond(Bond_ff(mol_ff.atoms[at1_idx],mol_ff.atoms[at2_idx],frc_val,length_val))

    return mol_ff

def save_parmed_angles(mol_ff,mol_pmd):

    for angle_pmd in mol_pmd.angles:
        at1_idx=angle_pmd.atom1.idx
        at2_idx=angle_pmd.atom2.idx
        at3_idx=angle_pmd.atom3.idx
        frc_val=angle_pmd.type.k
        ref_angle_val=angle_pmd.type.theteq

        if mol_ff.query_angle_atomIdx(at1_idx,at2_idx,at3_idx):
            angle_idx=mol_ff.query_angle_index(mol_ff.query_angle_atomIdx(at1_idx,at2_idx,at3_idx))
            mol_ff.angles[angle_idx].frc=frc_val
            mol_ff.angles[angle_idx].ref_angle=ref_angle_val

        else:
            mol_ff.add_angle(Angle_ff(mol_ff.atoms[at1_idx],mol_ff.atoms[at2_idx],mol_ff.atoms[at3_idx],frc_val,ref_angle_val))

    return mol_ff

def save_parmed_dihedrals(mol_ff,mol_pmd):

    for dihedral_pmd in mol_pmd.dihedrals:
        at1_idx=dihedral_pmd.atom1.idx
        at2_idx=dihedral_pmd.atom2.idx
        at3_idx=dihedral_pmd.atom3.idx
        at4_idx=dihedral_pmd.atom4.idx
        frc_val=float(dihedral_pmd.type.phi_k)
        period_val=int(dihedral_pmd.type.per)
        phase_val=round(float(dihedral_pmd.type.phase),1)
        improp=dihedral_pmd.improper

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

