# timber

import numpy as np
from .geometry import * 

##############################################################################

def setup_arrays(mol_ff):
    
    # return masks and array data for the mol_ff
    
    # bonds,frc_bond,length_bond,nonbond,a_ij,b_ij,e1_e2,q_ij,nb14,afep_angles,afep_dihedrals 
    
    # Init arrays

    # direct bonds
    bonds=np.zeros((len(mol_ff.atoms),len(mol_ff.atoms)),dtype=np.int8)
    
    # nonbonded mask for lj, elec
    nonbond=np.ones((len(mol_ff.atoms),len(mol_ff.atoms)),dtype=np.int8)
    np.fill_diagonal(nonbond, 0, wrap=False)

    # a_ij
    a_ij=np.zeros((len(mol_ff.atoms),len(mol_ff.atoms)),dtype=np.float32)
    b_ij=np.zeros((len(mol_ff.atoms),len(mol_ff.atoms)),dtype=np.float32)
    e1_e2=np.zeros((len(mol_ff.atoms),len(mol_ff.atoms)),dtype=np.float32)

    q_ij=np.zeros((len(mol_ff.atoms),len(mol_ff.atoms)),dtype=np.float32)

    frc_bond=np.zeros((len(mol_ff.atoms),len(mol_ff.atoms)),dtype=np.float32)
    length_bond=np.zeros((len(mol_ff.atoms),len(mol_ff.atoms)),dtype=np.float32)

    # nb14 mask
    nb14=np.zeros((len(mol_ff.atoms),len(mol_ff.atoms)),dtype=np.int8)

    # excluded 13
    remove13=np.zeros((len(mol_ff.atoms),len(mol_ff.atoms)),dtype=np.int8)

    # Init afep atom dicts
    afep_angles=[]
    
    afep_dihedrals=[]

    for i in range(0,len(mol_ff.atoms)):
        afep_angles.append(np.zeros(len(mol_ff.angles),dtype=np.int8))
        afep_dihedrals.append(np.zeros(len(mol_ff.dihedrals),dtype=np.int8))
    
    # fill arrays
    for i in range(0,len(mol_ff.atoms)):
        for j in range(0,len(mol_ff.atoms)):
            if i<j:

                r1=mol_ff.atoms[i].lj_R
                r2=mol_ff.atoms[j].lj_R

                e1=mol_ff.atoms[i].lj_E
                e2=mol_ff.atoms[j].lj_E

                q1=mol_ff.atoms[i].atom_charge*18.2223
                q2=mol_ff.atoms[j].atom_charge*18.2223

                a_ij[i][j]=np.sqrt(e1*e2)*np.power((r1+r2),12)
                a_ij[j][i]=np.sqrt(e1*e2)*np.power((r1+r2),12)

                b_ij[i][j]=2*np.sqrt(e1*e2)*np.power((r1+r2),6)
                b_ij[j][i]=2*np.sqrt(e1*e2)*np.power((r1+r2),6)

                e1_e2[i][j]=e1*e2
                e1_e2[j][i]=e1*e2

                q_ij[i][j]=q1*q2
                q_ij[j][i]=q1*q2

    # save bond atom indices
    for bond in mol_ff.bonds:
        at1=bond.atom1.idx
        at2=bond.atom2.idx

        bonds[at1][at2]=1
        bonds[at2][at1]=1

        frc_bond[at1][at2]=bond.frc
        frc_bond[at2][at1]=bond.frc

        length_bond[at1][at2]=bond.length
        length_bond[at2][at1]=bond.length
            
    # angles
    ctr=0
    for angle in mol_ff.angles:
        at1=angle.atom1.idx
        at2=angle.atom2.idx
        at3=angle.atom3.idx
        
        afep_angles[at1][ctr]=1
        afep_angles[at2][ctr]=1
        afep_angles[at3][ctr]=1

        remove13[at1][at3]=1
        remove13[at3][at1]=1

        ctr+=1

    # dihedrals
    ctr=0
    for dihed in mol_ff.dihedrals:
        at1=dihed.atom1.idx
        at2=dihed.atom2.idx
        at3=dihed.atom3.idx
        at4=dihed.atom4.idx

        afep_dihedrals[at1][ctr]=1
        afep_dihedrals[at2][ctr]=1
        afep_dihedrals[at3][ctr]=1
        afep_dihedrals[at4][ctr]=1

        # bonds,angles13 due to improper
        if remove13[at1][at4]==0 and bonds[at1][at4]==0:
            nb14[at1][at4]=1
            nb14[at4][at1]=1

        ctr+=1

    # turn off direct bonds in nonbond
    nonbond=nonbond*(1-bonds)

    # turn off 1-4 terms in nonbond
    nonbond=nonbond*(1-nb14)

    # turn off 1-3 terms in nonbond
    nonbond=nonbond*(1-remove13)
    
    del remove13
        
    return bonds,frc_bond,length_bond,nonbond,a_ij,b_ij,e1_e2,q_ij,nb14,afep_angles,afep_dihedrals

def calc_bond_array(bonds,frc_bond,length_bond,pair_dist):
       
    term1=np.subtract(pair_dist,length_bond)
    term2=np.multiply(term1,term1,out=term1)
    
    ene_bond=np.multiply(frc_bond,term2,out=term2)
    
    return np.nan_to_num(np.multiply(bonds,ene_bond,out=ene_bond))

def calc_lj_array(nonbond,a_ij,b_ij,pair_dist):
    
    term1=np.divide(a_ij,np.power(pair_dist,12))
    term2=np.divide(b_ij,np.power(pair_dist,6))
    
    eneLJ=np.subtract(term1,term2,out=term1)
        
    return np.nan_to_num(np.multiply(nonbond,eneLJ,out=eneLJ),posinf=0,neginf=0)
    
def calc_elec_array(nonbond,q_ij,pair_dist):
        
    eneElec=np.divide(q_ij,pair_dist)
        
    return np.nan_to_num(np.multiply(nonbond,eneElec,out=eneElec),posinf=0,neginf=0)
    
def get_angle_ene(mol_ff,xyz):

    # angles
    ene_angle=[]

    for angle in mol_ff.angles:
        at1=angle.atom1.idx
        at2=angle.atom2.idx
        at3=angle.atom3.idx

        deg=get_angle(xyz[at1],xyz[at2],xyz[at3])

        # if nan from mdtraj
        if math.isnan(deg):
            deg=angle.ref_angle

        ene=(angle.frc*(np.power((np.radians(deg)-np.radians(angle.ref_angle)),2)))

        ene_angle.append(ene)

    return ene_angle

def get_dihed_ene(mol_ff,xyz):

    # dihedrals
    ene_dihed=[]

    for dihed in mol_ff.dihedrals:
        at1=dihed.atom1.idx
        at2=dihed.atom2.idx
        at3=dihed.atom3.idx
        at4=dihed.atom4.idx

        deg=get_dihedral(xyz[at1],xyz[at2],xyz[at3],xyz[at4])

        # if nan
        if math.isnan(deg):
            deg=phase

        ene=0.0
        for i in range(0,dihed.n_terms):
            frc=dihed.frc[i]
            period=dihed.period[i]
            phase=dihed.phase[i]

            ene+=((frc)*(1+np.cos(period*np.radians(deg)-np.radians(phase))))
            
        ene_dihed.append(ene)
        
    return ene_dihed

def return_energy_decomposition(mol_ff):

    bonds,frc_bond,length_bond,nonbond,a_ij,b_ij,e1_e2,q_ij,nb14,afep_angles,afep_dihedrals=setup_arrays(mol_ff)

    coordinates=[atom.xyz for atom in mol_ff.atoms]

    pair_dist=[]
    for i in range(0,len(mol_ff.atoms)):
        for j in range(0,len(mol_ff.atoms)):
            dist=cart_distance(mol_ff.atoms[i].xyz,mol_ff.atoms[j].xyz)
            pair_dist.append(dist)

    pair_dist=np.nan_to_num(np.reshape(np.array(pair_dist),(len(mol_ff.atoms),len(mol_ff.atoms))))

    bonds=0.5*np.sum(calc_bond_array(bonds,frc_bond,length_bond,pair_dist))
    angles=np.sum(get_angle_ene(mol_ff,coordinates))
    dihedrals=np.sum(get_dihed_ene(mol_ff,coordinates))
    lj14=(1/2.0)*0.5*np.sum(calc_lj_array(nb14,a_ij,b_ij,pair_dist))
    elec14= (1/1.2)*0.5*np.sum(calc_elec_array(nb14,q_ij,pair_dist))
    lj=0.5*np.sum(calc_lj_array(nonbond,a_ij,b_ij,pair_dist))
    elec=0.5*np.sum(calc_elec_array(nonbond,q_ij,pair_dist))

    return bonds,angles,dihedrals,lj14,elec14,lj,elec

