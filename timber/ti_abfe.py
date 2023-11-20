# timber

import os
import glob
from rdkit import Chem
from .geometry import find_centroid,Coord,rank_closest,rank_further,cart_distance,get_angle,get_dihedral
from .molecule_ff import Molecule_ff,get_rdkit_info

##############################################################################

def write_abfe_disang(df,protocol,res_seq,dir_1_name='core'):

    # submit all dir lig1->lig2
    for index,row in df.iterrows():
        name1_col=list(df.columns)[0]
        name1=row[name1_col]
        pair_dir=name1

        os.chdir(pair_dir)

        LIG=Molecule_ff(name='LIG')
        mol=Chem.SDMolSupplier(dir_1_name+'/LIG.sdf',removeHs=False,sanitize=False)[0]
        LIG=get_rdkit_info(mol,LIG)

        lig_centroid=find_centroid([atom.xyz for atom in LIG.atoms])
        n_close,n_close_idx,ca_close,ca_close_idx,c_close,c_close_idx=coord_protein_backbone(lig_centroid,'complex/complex_ligands.pdb',res_seq+1)
        
        lig_atom1_idx,lig_atom2_idx,lig_atom3_idx,heavy_atom_position=closest_ligand_atoms(mol,lig_centroid)

        # Get the dist,angle,torsion values and write files
        dist1=cart_distance(heavy_atom_position[lig_atom1_idx],n_close)
        angle1=get_angle(heavy_atom_position[lig_atom1_idx],n_close,ca_close)
        angle2=get_angle(heavy_atom_position[lig_atom1_idx],heavy_atom_position[lig_atom2_idx],n_close)
        tor1=get_dihedral(c_close,ca_close,n_close,heavy_atom_position[lig_atom1_idx])
        tor2=get_dihedral(ca_close,n_close,heavy_atom_position[lig_atom1_idx],heavy_atom_position[lig_atom2_idx])
        tor3=get_dihedral(n_close,heavy_atom_position[lig_atom1_idx],heavy_atom_position[lig_atom2_idx],heavy_atom_position[lig_atom3_idx])

        write_disang('disang.RST',lig_atom1_idx+1,lig_atom2_idx+1,lig_atom3_idx+1,n_close_idx+int(len(mol.GetAtoms())),ca_close_idx+int(len(mol.GetAtoms())),c_close_idx+int(len(mol.GetAtoms())),dist1,angle1,angle2,tor1,tor2,tor3)

        write_disang('disang.shift.RST',lig_atom1_idx+1+int(len(mol.GetAtoms())),lig_atom2_idx+1+int(len(mol.GetAtoms())),lig_atom3_idx+1+int(len(mol.GetAtoms())),n_close_idx+2*int(len(mol.GetAtoms())),ca_close_idx+2*int(len(mol.GetAtoms())),c_close_idx+2*int(len(mol.GetAtoms())),dist1,angle1,angle2,tor1,tor2,tor3)

        if protocol=='absolute-three-step':
            write_disang('disang.pair.RST',lig_atom1_idx+1,lig_atom2_idx+1,lig_atom3_idx+1,n_close_idx+2*int(len(mol.GetAtoms())),ca_close_idx+2*int(len(mol.GetAtoms())),c_close_idx+2*int(len(mol.GetAtoms())),dist1,angle1,angle2,tor1,tor2,tor3,pair_offset=int(len(mol.GetAtoms())))

        os.chdir('../')

def coord_protein_backbone(lig_centroid,pdb_file,res_seq):

    # get all protein carbon alpha positions
    prot_calpha={}
    prot_points=[]
    with open(pdb_file,'r') as f:
        for line in f:
             if (len(line.split())>1 and (line.split()[0]=='ATOM') and (line.split()[2]=='CA')):
                 x=float(line.split()[5])
                 y=float(line.split()[6])
                 z=float(line.split()[7])

                 resi=int(line.split()[4])
                 prot_calpha.update({resi:Coord(x,y,z)})
                 prot_points.append(Coord(x,y,z))

    # amber format residue number to use as the anchor
    resi_closest=int(res_seq)

    # Get the coords of CA, N, C for closest residue
    ca_close=None
    ca_close_idx=None
    n_close=None
    n_close_idx=None
    c_close=None
    c_close_idx=None
    with open(pdb_file,'r') as f:
        for line in f:
            if line.split()[0]=='ATOM':
                if int(line.split()[4])==resi_closest:
                    x=float(line.split()[5])
                    y=float(line.split()[6])
                    z=float(line.split()[7])
                    if line.split()[2]=='CA':
                        ca_close=Coord(x,y,z)
                        ca_close_idx=int(line.split()[1])
                    elif line.split()[2]=='C':
                        c_close=Coord(x,y,z)
                        c_close_idx=int(line.split()[1])
                    elif line.split()[2]=='N':
                        n_close=Coord(x,y,z)
                        n_close_idx=int(line.split()[1])

    return n_close,n_close_idx,ca_close,ca_close_idx,c_close,c_close_idx

def closest_ligand_atoms(rdmol,lig_centroid):

    # Ligand atom closest to centroid
    heavy_atom_position=[]
    for i in range(0,len(rdmol.GetAtoms())):
        if rdmol.GetAtomWithIdx(i).GetAtomicNum()!=1:
            pos=rdmol.GetConformer().GetAtomPosition(i)
            heavy_atom_position.append(Coord(float(pos.x),float(pos.y),float(pos.z)))

    lig_atom1_idx=rank_closest(lig_centroid,heavy_atom_position)
    lig_atom2_idx=rank_further(heavy_atom_position[lig_atom1_idx],heavy_atom_position)

    lig_angles=[]
    for i in range(0,len(heavy_atom_position)):
        if (i!=lig_atom1_idx) and (i!=lig_atom2_idx):
            my_angle=get_angle(heavy_atom_position[i],heavy_atom_position[lig_atom1_idx],heavy_atom_position[lig_atom2_idx])
            if my_angle<150:
                lig_angles.append(my_angle)
            else:
                lig_angles.append(0)
        else:
            lig_angles.append(0)

    lig_atom3_idx=int(lig_angles.index(max(lig_angles)))

    return lig_atom1_idx,lig_atom2_idx,lig_atom3_idx,heavy_atom_position

def write_disang(disang_file_name,lig_atom1_idx,lig_atom2_idx,lig_atom3_idx,n_close_idx,ca_close_idx,c_close_idx,dist1,angle1,angle2,tor1,tor2,tor3,pair_offset=None):

    with open(disang_file_name,'w') as f:
        f.write('# Harmonic restraints\n')
        f.write('&rst iat=%s,%d, r1=0., r2=%6.2f, r3=%6.2f, r4=%6.2f, rk2=1000.0, rk3=1000.0,/\n' % (lig_atom1_idx,n_close_idx,dist1,dist1,dist1+10.0))
        if pair_offset:
            f.write('&rst iat=%s,%d, r1=0., r2=%6.2f, r3=%6.2f, r4=%6.2f, rk2=1000.0, rk3=1000.0,/\n' % (lig_atom1_idx+pair_offset,n_close_idx,dist1,dist1,dist1+10.0))

        if angle1>0:
            f.write('&rst iat=%d,%d,%d, r1=0., r2=%6.2f, r3=%6.2f, r4=180., rk2=10.0, rk3=10.0,/\n'  % (lig_atom1_idx,n_close_idx,ca_close_idx,angle1,angle1))
            if pair_offset:
                f.write('&rst iat=%d,%d,%d, r1=0., r2=%6.2f, r3=%6.2f, r4=180., rk2=10.0, rk3=10.0,/\n'  % (lig_atom1_idx+pair_offset,n_close_idx,ca_close_idx,angle1,angle1))
        else:
            f.write('&rst iat=%d,%d,%d, r1=-180., r2=%6.2f, r3=%6.2f, r4=0., rk2=10.0, rk3=10.0,/\n'  % (lig_atom1_idx,n_close_idx,ca_close_idx,angle1,angle1))
            if pair_offset:
                f.write('&rst iat=%d,%d,%d, r1=-180., r2=%6.2f, r3=%6.2f, r4=0., rk2=10.0, rk3=10.0,/\n'  % (lig_atom1_idx+pair_offset,n_close_idx,ca_close_idx,angle1,angle1))

        if angle2>0:
            f.write('&rst iat=%d,%d,%d, r1=0., r2=%6.2f, r3=%6.2f, r4=180., rk2=10.0, rk3=10.0,/\n'  % (lig_atom1_idx,lig_atom2_idx,n_close_idx,angle2,angle2))
            if pair_offset:
                f.write('&rst iat=%d,%d,%d, r1=0., r2=%6.2f, r3=%6.2f, r4=180., rk2=10.0, rk3=10.0,/\n'  % (lig_atom1_idx+pair_offset,lig_atom2_idx+pair_offset,n_close_idx,angle2,angle2))
        else:
            f.write('&rst iat=%d,%d,%d, r1=-180., r2=%6.2f, r3=%6.2f, r4=0., rk2=10.0, rk3=10.0,/\n'  % (lig_atom1_idx,lig_atom2_idx,n_close_idx,angle2,angle2))
            if pair_offset:
                f.write('&rst iat=%d,%d,%d, r1=-180., r2=%6.2f, r3=%6.2f, r4=0., rk2=10.0, rk3=10.0,/\n'  % (lig_atom1_idx+pair_offset,lig_atom2_idx+pair_offset,n_close_idx,angle2,angle2))

        if tor1>0:
            f.write('&rst iat=%d,%d,%d,%d, r1=0., r2=%6.2f, r3=%6.2f, r4=180., rk2=10.0, rk3=10.0,/\n' % (c_close_idx,ca_close_idx,n_close_idx,lig_atom1_idx,tor1,tor1))
            if pair_offset:
                f.write('&rst iat=%d,%d,%d,%d, r1=0., r2=%6.2f, r3=%6.2f, r4=180., rk2=10.0, rk3=10.0,/\n' % (c_close_idx,ca_close_idx,n_close_idx,lig_atom1_idx+pair_offset,tor1,tor1))
        else:
            f.write('&rst iat=%d,%d,%d,%d, r1=-180., r2=%6.2f, r3=%6.2f, r4=0., rk2=10.0, rk3=10.0,/\n' % (c_close_idx,ca_close_idx,n_close_idx,lig_atom1_idx,tor1,tor1))
            if pair_offset:
                f.write('&rst iat=%d,%d,%d,%d, r1=-180., r2=%6.2f, r3=%6.2f, r4=0., rk2=10.0, rk3=10.0,/\n' % (c_close_idx,ca_close_idx,n_close_idx,lig_atom1_idx+pair_offset,tor1,tor1))

        if tor2>0:
            f.write('&rst iat=%d,%d,%d,%d, r1=0., r2=%6.2f, r3=%6.2f, r4=180., rk2=10.0, rk3=10.0,/\n' % (ca_close_idx,n_close_idx,lig_atom1_idx,lig_atom2_idx,tor2,tor2))
            if pair_offset:
                f.write('&rst iat=%d,%d,%d,%d, r1=0., r2=%6.2f, r3=%6.2f, r4=180., rk2=10.0, rk3=10.0,/\n' % (ca_close_idx,n_close_idx,lig_atom1_idx+pair_offset,lig_atom2_idx+pair_offset,tor2,tor2))
        else:
            f.write('&rst iat=%d,%d,%d,%d, r1=-180., r2=%6.2f, r3=%6.2f, r4=0., rk2=10.0, rk3=10.0,/\n' % (ca_close_idx,n_close_idx,lig_atom1_idx,lig_atom2_idx,tor2,tor2))
            if pair_offset:
                f.write('&rst iat=%d,%d,%d,%d, r1=-180., r2=%6.2f, r3=%6.2f, r4=0., rk2=10.0, rk3=10.0,/\n' % (ca_close_idx,n_close_idx,lig_atom1_idx+pair_offset,lig_atom2_idx+pair_offset,tor2,tor2))

        if tor3>0:
            f.write('&rst iat=%d,%d,%d,%d, r1=0., r2=%6.2f, r3=%6.2f, r4=180., rk2=10.0, rk3=10.0,/\n' % (n_close_idx,lig_atom1_idx,lig_atom2_idx,lig_atom3_idx,tor3,tor3))
            if pair_offset:
                f.write('&rst iat=%d,%d,%d,%d, r1=0., r2=%6.2f, r3=%6.2f, r4=180., rk2=10.0, rk3=10.0,/\n' % (n_close_idx,lig_atom1_idx+pair_offset,lig_atom2_idx+pair_offset,lig_atom3_idx+pair_offset,tor3,tor3))
        else:
            f.write('&rst iat=%d,%d,%d,%d, r1=-180., r2=%6.2f, r3=%6.2f, r4=0., rk2=10.0, rk3=10.0,/\n' % (n_close_idx,lig_atom1_idx,lig_atom2_idx,lig_atom3_idx,tor3,tor3))
            if pair_offset:
                f.write('&rst iat=%d,%d,%d,%d, r1=-180., r2=%6.2f, r3=%6.2f, r4=0., rk2=10.0, rk3=10.0,/\n' % (n_close_idx,lig_atom1_idx+pair_offset,lig_atom2_idx+pair_offset,lig_atom3_idx+pair_offset,tor3,tor3))

