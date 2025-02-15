#!/usr/bin/env python
import os
import sys
import argparse
import mdtraj as md
from boresch_restraints import *

## Command line arguments ##
parser = argparse.ArgumentParser(description='Return protein and ligand atom indices for Boresch restraints (zero index)\n')
parser.add_argument('-i',help='trajectory file',dest='i',default='image.nc',required=False)
parser.add_argument('-p',help='prmtop file',dest='p',default='prot_UNL.prmtop',required=False)
parser.add_argument('-lig',help='ligand mask string',dest='lig',default='UNL',required=False)
parser.add_argument('-sdf',help='ligand SDF file',dest='sdf',default='UNL.sdf',required=False)

args=vars(parser.parse_args())

## Settings ##
prmtopfile=args['p']
traj_file=args['i']
ligand_file=args['sdf']
ligand_str=args['lig']

## Load traj ##
traj_md = md.load(traj_file,top=prmtopfile,stride=1)

## Get 3 ligand atoms ##
l1,l2,l3,lig_len=select_ligand_atoms(lig=ligand_file, traj=traj_md, ligand=ligand_str)

## Get 3 protein atoms ##
p1,p2,p3,l1,l2,l3=select_Boresch_atoms(traj=traj_md, mol2_lig=ligand_file, ligand_atoms = [l1,l2,l3], ligand=ligand_str)[0]

lig_start=traj_md.topology.select('resn %s' % (ligand_str))[0]
prot_start=traj_md.topology.select('protein')[0]

## Reset the atom indices ##
p1=p1-prot_start
p2=p2-prot_start
p3=p3-prot_start

## Reset the ligand indices ##
l1=l1-lig_start
l2=l2-lig_start
l3=l3-lig_start

## Print result ##
print('Protein: %d,%d,%d' % (p1,p2,p3))
print('Ligand: %d,%d,%d' % (l1,l2,l3))

