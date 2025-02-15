#!/usr/bin/env python
import argparse
import sys
import os
import numpy as np
import time

from ase_interface import ANIENS
from ase_interface import aniensloader

import ase
from ase import units
from ase.io import read, write
from ase.optimize import BFGS, LBFGS
from ase.constraints import FixInternals

from rdkit import Chem
from rdkit.Chem.PropertyMol import *
from rdkit.Chem import SDWriter

###############################################################################
def load_ani(mol_xyz,gpu_id):
	mol = read(mol_xyz)
	
	# Set ANI calculator
	mol.set_calculator(ANIENS(aniensloader('/usr/prog/cadd/amber_tools/alchemistry/ASE_ANI/ani_models/ani-2x_8x.info',gpu_id)))

	# check for accepted atom types
	accepted_types=['C','H','N','O','S','Cl','F']
	for i,atom in enumerate(mol):
		if str(atom.symbol) not in accepted_types:
			raise ValueError('Only HCNOSFCl atom types permitted!')

	return mol

def run_single_point(mol_xyz,gpu_id):
	mol=load_ani(mol_xyz,gpu_id)

	e=mol.get_potential_energy()

	# atomic energy breakdown
	atomic_breakdown=np.mean(mol.calc.get_atomicenergies(sae=False),axis=0)

	return e,atomic_breakdown 

def run_optimize(mol_xyz,traj_out_file,gpu_id):
	mol=load_ani(mol_xyz,gpu_id)

	dyn = LBFGS(mol)
	dyn.run(fmax=0.001)

	write(traj_out_file, mol, format='xyz', append=True)
	e=mol.get_potential_energy()

	atomic_breakdown=np.mean(mol.calc.get_atomicenergies(sae=False),axis=0)

	return e,atomic_breakdown

def run_torsion_scan(rd_mol,mol_xyz,torsion_list,n_steps,step_size,gpu_id):
	# ERROR CHECKS
	assert isinstance(torsion_list,list), 'Provide list of atom indices!'

	# ASE expects zero indexed values
	torsion_list=[int(i)-1 for i in torsion_list]

	if len(torsion_list)!=4:
		raise ValueError('Must provide 4 atom indices!')

	mol=load_ani(mol_xyz,gpu_id)

	ene_values=[]
	dihed_values=[]
	sd_mols=[]

	# RUN CODE
	for i in range(0,n_steps):
		if i==0:
			mol.rotate_dihedral(a1=torsion_list[0],a2=torsion_list[1],a3=torsion_list[2],a4=torsion_list[3], angle=float(0))
		else:
			mol.rotate_dihedral(a1=torsion_list[0],a2=torsion_list[1],a3=torsion_list[2],a4=torsion_list[3], angle=float(step_size))
		dihedral1 = [mol.get_dihedral(*torsion_list) * 3.141592653589793 / 180,torsion_list]
		c = FixInternals(bonds=[], angles=[],dihedrals=[dihedral1])
		mol.set_constraint(c)
		dyn = LBFGS(mol)
		dyn.run(fmax=0.001)

		ene_values.append(mol.get_potential_energy())
		dihed_values.append(dihedral1[0]* (180/3.141592653589793))

		write('opt_traj.xyz', mol, format='xyz', append=False)
		my_rd_mol=Chem.Mol(rd_mol)
		my_rd_mol=update_cartesian(my_rd_mol,'opt_traj.xyz')
		os.system('rm opt_traj.xyz')

		sd_mols.append(my_rd_mol)

	return dihed_values,ene_values,sd_mols

###############################################################################

def mol_type(mol_input):
	if not os.path.isfile(mol_input):
		print('\nCannot find input!\n')
		sys.exit()

	accepted=['sdf','mol2']
	mol_type=mol_input.split('.')[-1]
	if mol_type not in accepted:
		print('\nUnrecognized mol input format!\n')
		sys.exit()

	return mol_type

def convert_sdf(rdmol):
	with open('tmpAni.xyz','w') as f_out:
		f_out.write('%d\n' % (len(rdmol.GetAtoms())))
		f_out.write('comment\n')
		for i in range(0,len(rdmol.GetAtoms())):
			pos=rdmol.GetConformer().GetAtomPosition(i)
			f_out.write('%s %lf %lf %lf\n' % (rdmol.GetAtomWithIdx(i).GetSymbol(),pos.x,pos.y,pos.z))

def rd_write(sd_out,mol_out):
	# Sort mol list by eV
	mol_out.sort(key=lambda m: float(m.GetProp('ANI_eV')), reverse=False)
	base_ene=float(mol_out[0].GetProp('ANI_eV'))

	# Set the kcal values
	for mol in mol_out:
		val=(float(mol.GetProp('ANI_eV'))-base_ene)*23.061
		mol.SetProp('ANI_kcal',val)

	# Write SDF file
	writer=SDWriter(sd_out)
	for mol in mol_out:
		writer.write(mol)
	writer.flush()

def update_cartesian(rd_mol,xyz_mol):
	with open(xyz_mol,'r') as f:
		f.readline()
		f.readline()
		data=f.readlines()
		for i in range(0,len(rd_mol.GetAtoms())):
			x=float(data[i].split()[1])
			y=float(data[i].split()[2])
			z=float(data[i].split()[3])

			rd_mol.GetConformer().SetAtomPosition(i,(x,y,z))
	return rd_mol 

def parse_torsion_output(sdf_input,dihed_values,ene_values,rd_scan_mols,torsion_list,name=None):

	torsion_list=[int(i) for i in torsion_list]

	for i in range(0,len(dihed_values)):
		rd_scan_mols[i].SetProp('Dihedral angle',str(dihed_values[i]))	
		rd_scan_mols[i].SetProp('ANI_eV',str(ene_values[i]))

		val=(ene_values[i]-min(ene_values))*23.061
		rd_scan_mols[i].SetProp('ANI_kcal',str(val))

	if name==None:
		sd_out=args['i'].split('.')[0]+'_ani-dihedral.sdf'
		data_out=args['i'].split('.')[0]+'_ani-dihedral.dat'
	else:
		sd_out=name+'_ani-dihedral.sdf'
		data_out=name+'_ani-dihedral.dat'

	with open(data_out,'w') as f_out:
		for i in range(0,len(dihed_values)):
			f_out.write('%lf %lf\n' % (dihed_values[i],(ene_values[i]-min(ene_values))*23.061))
	
	# Write SDF file
	writer=SDWriter(sd_out)
	for mol in rd_scan_mols:
		writer.write(mol)
	writer.flush()

	# Sort the mols by dihedral and write out sort file
	rd_scan_mols.sort(key=lambda m: float(m.GetProp('Dihedral angle')), reverse=False)

	if name==None:
		sd_out=args['i'].split('.')[0]+'_ani-dihedral_sort.sdf'
		data_out=args['i'].split('.')[0]+'_ani-dihedral_sort.dat'
	else:
		sd_out=name+'_ani-dihedral_sort.sdf'
		data_out=name+'_ani-dihedral_sort.dat'

	with open(data_out,'w') as f_out:
		for i in range(0,len(dihed_values)):
			f_out.write('%lf %lf\n' % (float(rd_scan_mols[i].GetProp('Dihedral angle')),float(rd_scan_mols[i].GetProp('ANI_kcal'))))

	# Write SDF file
	writer=SDWriter(sd_out)
	for mol in rd_scan_mols:
		writer.write(mol)
	writer.flush()

	# Write a pml file to load into pymol
	with open('load_dihe.pml','w') as f_out:
		f_out.write('load %s\n' % (sd_out))

		for i in range(0,len(dihed_values)):
			f_out.write('select val, (%s and id %d and state %d)\n' % (sd_out.split('.')[0],torsion_list[1],i))
			f_out.write('label val, "%5.2f deg %5.2f kcal/mol" \n' % (float(rd_scan_mols[i].GetProp('Dihedral angle')),float(rd_scan_mols[i].GetProp('ANI_kcal'))))

		f_out.write('set label_position, (2,2,2)\n')
		f_out.write('delete val\n')
		f_out.write('select dihe, (%s and id %d,%d,%d,%d)\n' % (sd_out.split('.')[0],torsion_list[0],torsion_list[1],torsion_list[2],torsion_list[3]))
		f_out.write('show spheres, dihe\n')
		f_out.write('set sphere_scale, 0.2\n')

def perform_multi(args,rd_suppl,gpu_id,name=None):
	# Single-point energy
	if args['o']==False and args['s']==None:
		mol_out=[]

		for mol in rd_suppl:
			pm=PropertyMol(mol)
			
			convert_sdf(mol)
			e,b=run_single_point('tmpAni.xyz',gpu_id)
			os.system('rm tmpAni.xyz')

			pm.SetProp('ANI_eV',e)
			mol_out.append(pm)

			if args['break']:
				for i in range(0,len(pm.GetAtoms())):
					print('%d %lf' % (int(pm.GetAtomWithIdx(i).GetAtomicNum()),b[i]))

			print('Single-point (eV): %lf' % (e))

		if name==None:
			sd_out_name=args['i'].split('.')[0]+'_ani-single.sdf'
		else:
			sd_out_name=name+'_ani-single.sdf'
		rd_write(sd_out_name,mol_out)

	# Optimization
	if args['o']==True and args['s']==None:
		mol_out=[]

		for mol in rd_suppl:
			pm=PropertyMol(mol)

			convert_sdf(mol)
			e,b=run_optimize('tmpAni.xyz','opt_traj.xyz',gpu_id)

			pm=update_cartesian(pm,'opt_traj.xyz')
			os.system('rm opt_traj.xyz tmpAni.xyz')

			pm.SetProp('ANI_eV',e)
			mol_out.append(pm)

			if args['break']:
				for i in range(0,len(pm.GetAtoms())):
					print('%d %lf' % (int(pm.GetAtomWithIdx(i).GetAtomicNum()),b[i]))

			print('Optimized (eV): %lf' % (e))

		if name==None:
			sd_out_name=args['i'].split('.')[0]+'_ani-opt.sdf'
		else:
			sd_out_name=name+'_ani-opt.sdf'
		rd_write(sd_out_name,mol_out)

	# Torsion scan
	if args['o']==False and args['s']!=None:
		if len(rd_suppl)==1:
			print('running dihedral scan %d %d %d %d ... \n' % (int(args['s'][0]),int(args['s'][1]),int(args['s'][2]),int(args['s'][3])))

			pm=PropertyMol(Chem.Mol(rd_suppl[0]))

			convert_sdf(pm)
			dihed_values,ene_values,rd_scan_mols=run_torsion_scan(pm,'tmpAni.xyz',args['s'],37,10,gpu_id)
			parse_torsion_output(args['i'],dihed_values,ene_values,rd_scan_mols,args['s'],name)

			os.system('rm tmpAni.xyz')

		else:
			print('SDF file contains >1 mol!')
			sys.exit()

	if args['o']==True and args['s']!=None:
		print('Cannot opt and scan')
		sys.exit()

###############################################################################

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='ANI model to optimize/scan sdf file. NOTE: torsion indices as pymol->ID\n')
	parser.add_argument('-i',help='Input mol file',required=True)
	parser.add_argument('-o',help='Optimize flag',action='store_true')
	parser.add_argument('-n',help='Output name',required=False)
	parser.add_argument('-break',help='Breakdown of atomic energies',action='store_true')
	parser.add_argument('-s',help='Scan dihedral flag',action='store',nargs=4)
	parser.add_argument('-gpu',help='GPU ID',type=int,required=False,default=0)
	args=vars(parser.parse_args())

	gpu_id=int(args['gpu'])

	if mol_type(args['i'])=='sdf':
		suppl=Chem.SDMolSupplier(args['i'],removeHs=False)

		perform_multi(args,suppl,gpu_id,args['n'])

	elif mol_type(args['i'])=='mol2':
		# hack way to deal with mol2
		suppl=[]
		mol=Chem.MolFromMol2File(args['i'],removeHs=False)
		suppl.append(mol)

		perform_multi(args,suppl,gpu_id,args['n'])
