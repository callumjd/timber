#!/usr/bin/env python
import argparse
import glob
import os
import sys
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem,EditableMol
from rdkit.Chem import rdMolTransforms,rdmolops
from rdkit.Chem import SDWriter

#
# Callum Dickson, NIBR CADD: callum.dickson@novartis.com
#

###############################################################################
## Routines ##
###############################################################################

def check_file(file_in):
    output=False
    if os.path.exists(file_in) and os.path.getsize(file_in)>0:
        output=True
    return output

class Atom_rocs(object):
    def __init__(self,atomic_num,x=None,y=None,z=None,dist=None):
        self.atomic_num=atomic_num
        self.x=x
        self.y=y
        self.z=z
        self.dist=dist

# x,y,z cartesian coordinate object
class Coord(object):
    def __init__(self,x,y,z):
        self.x=x
        self.y=y
        self.z=z

## distance ##
def cart_distance(coord1,coord2):
    x1=coord2.x-coord1.x
    y1=coord2.y-coord1.y
    z1=coord2.z-coord1.z

    tot=(x1*x1)+(y1*y1)+(z1*z1)
    return math.sqrt(tot)

def parse_watmd(input_file,acceptor_cutoff,donor_cutoff,occupancy_cutoff):

    # returns a list of atoms with x,y,z coordinates
    atom_list=[]
    with open(input_file,'r') as f:
        for line in f:
            if len(line.split())>1:
                if (line.split()[0]=='ATOM' or line.split()[0]=='HETATM'):
                    if str(line[17:21].strip())=='CVX':
                        x=float(line[30:38])
                        y=float(line[38:46])
                        z=float(line[46:54])
                        chg=float(line[55:62])
                        occ=float(line[62:69])

                        # make a rocs atom
                        if occ>occupancy_cutoff:
                            # acceptor Ar
                            if chg<acceptor_cutoff:
                                atom=Atom_rocs(18,x,y,z)
                            # donor Kr
                            elif chg>donor_cutoff:
                                atom=Atom_rocs(36,x,y,z)
                            # grease Xe
                            else:
                                atom=Atom_rocs(54,x,y,z)

                            atom_list.append(atom)
    return atom_list

def parse_site(prot,lig,site_cutoff):

    # use ligand SDF file to define the binding site excluded volume
    if lig.split('.')[-1]=='sdf':
        mol=Chem.SDMolSupplier(lig)[0]
        if not mol:
            print('Error: could not load ligand SDF!\n')
            sys.exit()

        ligand_atoms=[]
        for i in range(0,len(mol.GetAtoms())):
            pos=mol.GetConformer().GetAtomPosition(i)
            ligand_atoms.append(Coord(pos.x,pos.y,pos.z))

    # use the water query to define the binding site excluded volume
    elif lig.split('.')[-1]=='pdb':

        ligand_atoms=[]
        with open(lig,'r') as f:
            for line in f:
                if len(line.split())>1:
                    # save all water sites 
                    if line.split()[0]=='ATOM':
                        x=float(line[30:38])
                        y=float(line[38:46])
                        z=float(line[46:54])
                        coord=Coord(x,y,z)

                        ligand_atoms.append(coord)

    else:
        print('Error: unrecognized ligand format!\n')
        sys.exit()

    protein_atoms=[]
    with open(prot,'r') as f:
        for line in f:
            if len(line.split())>1:
                # find CA atoms near any atom of the ligand
                if line.split()[0]=='ATOM' and line.split()[2] in ['CA']:
                    x=float(line[30:38])
                    y=float(line[38:46])
                    z=float(line[46:54])
                    coord=Coord(x,y,z)

                    for atom in ligand_atoms:
                        dist=cart_distance(atom,coord)
                        if dist<site_cutoff:
                            # Excluded atoms radon
                            exclude=Atom_rocs(86,x,y,z,dist)
                            protein_atoms.append(exclude)
                            break

    # sort by distance
    protein_atoms=sorted(protein_atoms, key=lambda x: float(x.dist), reverse=False)

    return protein_atoms

def create_rocs_mol(query_atoms):

    em=Chem.EditableMol(Chem.Mol())
    for at in query_atoms:
        em.AddAtom(Chem.Atom(at.atomic_num))

    mol=em.GetMol()
    conf=Chem.Conformer(len(query_atoms))
    mol.AddConformer(conf)

    for i in range(0,len(query_atoms)):
        x=query_atoms[i].x
        y=query_atoms[i].y
        z=query_atoms[i].z

        mol.GetConformer().SetAtomPosition(i,(x,y,z))

    return mol

def run(input_file,prot,lig,output_name): 

    # Cut-offs to classify water as acceptor, donor, or grease
    acceptor_cutoff=-0.5
    donor_cutoff=0.5
    occupancy_cutoff=0
    # grease will lie between the acceptor_cutoff and donor_cutoff

    # Halogens:
    # acceptor Ar
    # donor Kr
    # grease Xe
    # excluded volume Rn

    if not check_file(input_file):
        print('Error: cannot find %s!\n' % (input_file))
        sys.exit()
    if lig and not prot:
        print('Error: must specify both protein and ligand for excluded volume!\n')
        sys.exit()

    if prot:
        if not check_file(prot):
            print('Error: cannot find %s!\n' % (prot))
            sys.exit()
    if lig:
        if not check_file(lig):
            print('Error: cannot find %s!\n' % (lig))
            sys.exit()

    print('\nProcessing %s\n' % (input_file))

    # parse the watmd input
    query_atoms=parse_watmd(input_file,acceptor_cutoff,donor_cutoff,occupancy_cutoff)

    # parse the protein
    # WARNING: adding too many excluded volume points appears 
    # to give poor ROCS output
    # Apply hard cut-off of 20
    if prot and lig:
        # cut-off of 7 to carve the site
        print('Using %s and %s for binding site excluded volume\n' % (prot,lig))
        excluded_atoms=parse_site(prot,lig,7)
        query_atoms=query_atoms+excluded_atoms[0:20]
    elif prot and not lig:
        # cut-off of 7 to carve the site
        print('Using %s and %s for binding site excluded volume\n' % (prot,input_file))
        excluded_atoms=parse_site(prot,input_file,7)
        query_atoms=query_atoms+excluded_atoms[0:20]

    # make the output mol
    output_mol=create_rocs_mol(query_atoms)
    output_mol.SetProp('_Name',output_name)

    # save to output
    print('Writing %s\n' % (output_name))
    writer=SDWriter(output_name)
    writer.write(output_mol)
    writer.flush()

###############################################################################
## Parse arguments ##
###############################################################################

def main(argv=None):
    ## Command line arguments ##
    parser = argparse.ArgumentParser(description='Convert WATMD to ROCS query\n')

    parser.add_argument('-i','--input',dest='i',help='WATMD waters',type=str,required=True)
    parser.add_argument('-p','--prot',dest='p',help='Protein PDB file for excluded volume',type=str,required=False)
    parser.add_argument('-l','--ligand',dest='l',help='Ligand SDF file for excluded volume (if not specified, water query used)',type=str,required=False)
    parser.add_argument('-o','--output',help='Output name for sdf file',dest='o',default='water_query.sdf',required=False)

    args=vars(parser.parse_args())

###############################################################################
## Run selected protocols ##
###############################################################################
    run(input_file=args['i'],prot=args['p'],lig=args['l'],output_name=args['o'])

# MAIN
if __name__ == '__main__':
    main()

