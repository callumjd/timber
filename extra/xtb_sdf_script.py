#!/usr/bin/env python

#
# XTB Wrapper
# Callum Dickson, 8 May 2019
#
# Updated 12 June 2020 for XTB v6.3
#

import argparse
import sys
import os
import numpy as np
import time

from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem.PropertyMol import *
from rdkit.Chem import SDWriter
from rdkit.Chem import rdmolops

import distutils.spawn

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

def xtb_opt(rd_mol,mol_xyz,xtb_exe,gbsa_flag,formal_charge,output_name,freeze_list,freeze_file):

	if (freeze_list!=None or freeze_file!=None):
		with open('cntrl_xtb.in','w') as f:

			if (freeze_list!=None and freeze_file==None):
				freeze_list=[int(i) for i in freeze_list]
				f.write('$constrain\n')
				f.write('force constant=5\n')
				f.write('dihedral: %d,%d,%d,%d,auto\n' % (freeze_list[0],freeze_list[1],freeze_list[2],freeze_list[3]))


			elif (freeze_list==None and freeze_file!=None):
				with open(freeze_file,'r') as fl_file:
					for line in fl_file:
						if len(line.split('-'))>1:
							f.write('$constrain\n')
							f.write('force constant=5\n')
						
							fl_1=int(line.split('-')[0])
							fl_2=int(line.split('-')[1])
							fl_3=int(line.split('-')[2])
							fl_4=int(line.split('-')[3])

							f.write('dihedral: %d,%d,%d,%d,auto\n' % (fl_1,fl_2,fl_3,fl_4))
				
						elif len(line.split('-'))==1:
							f.write('$fix\n')
							at_1=int(line.strip())
							f.write('atoms: %d\n' % (at_1))
			f.write('$end\n')

	# run an optimization on a single sdf mol
	if gbsa_flag==None:
		if (freeze_list!=None or freeze_file!=None):
			os.system('%s %s --opt -c %d -I cntrl_xtb.in > xtb_raw_out' % (xtb_exe,mol_xyz,formal_charge))
			os.system('rm cntrl_xtb.in')
		else:
			os.system('%s %s --opt -c %d > xtb_raw_out' % (xtb_exe,mol_xyz,formal_charge))
	else:
		if (freeze_list!=None or freeze_file!=None):
			os.system('%s %s --opt -c %d --alpb %s -I cntrl_xtb.in > xtb_raw_out' % (xtb_exe,mol_xyz,formal_charge,gbsa_flag))
			os.system('rm cntrl_xtb.in')
		else:
			os.system('%s %s --opt -c %d --alpb %s > xtb_raw_out' % (xtb_exe,mol_xyz,formal_charge,gbsa_flag))

	ene_values=[]
	sd_mols=[]

	with open('xtbopt.xyz','r') as f:
		f.readline()
		ene_values.append(float(f.readline().split()[1].strip()))

       # convert xyz structures to sdf
	local_mol=Chem.Mol(rd_mol)

	with open('xtbopt.xyz','r') as f:
		f.readline()
		f.readline()

		counter=0
		for line in f:
			x=float(line.split()[1])
			y=float(line.split()[2])
			z=float(line.split()[3])

			local_mol.GetConformer().SetAtomPosition(counter,(x,y,z))
			counter+=1

	sd_mols.append(local_mol)

	if output_name is not None:
		output_name=str(output_name)
	else:
		output_name='opt_mols'


	writer=SDWriter(output_name+'_xtb.sdf')
	counter=0
	for mol in sd_mols:
		pm=PropertyMol(mol)
		pm.SetProp('xtb-hartree',ene_values[counter])
		writer.write(pm)
		counter+=1
	writer.flush()

	os.system('rm %s .xtboptok charges wbo xtbopt.log xtbopt.xyz xtbrestart xtb_raw_out xtbtopo.mol' % (mol_xyz))

def xtb_torsion_scan(rd_mol,mol_xyz,torsion_list,freeze_list,freeze_file,n_steps,step_size,xtb_exe,gbsa_flag,formal_charge,output_name):
	# ERROR CHECKS
	assert isinstance(torsion_list,list), 'Provide list of atom indices!'

	torsion_list=[int(i) for i in torsion_list]

	if len(torsion_list)!=4:
		raise ValueError('Must provide 4 atom indices!')

	if (freeze_list!=None and freeze_file==None):
		freeze_list=[int(i) for i in freeze_list]

	ene_values=[]
	dihed_values=[]
	sd_mols=[]

	start_angle=rdMolTransforms.GetDihedralDeg(rd_mol.GetConformer(),torsion_list[0]-1,torsion_list[1]-1,torsion_list[2]-1,torsion_list[3]-1)

	with open('cntrl_xtb.in','w') as f:
		f.write('$constrain\n')
		f.write('force constant=5\n')	
		f.write('dihedral: %d,%d,%d,%d,auto\n' % (torsion_list[0],torsion_list[1],torsion_list[2],torsion_list[3]))

		if (freeze_list!=None and freeze_file==None):
			f.write('dihedral: %d,%d,%d,%d,auto\n' % (freeze_list[0],freeze_list[1],freeze_list[2],freeze_list[3]))

		elif (freeze_list==None and freeze_file!=None):
			with open(freeze_file,'r') as fl_file:
				for line in fl_file:
					fl_1=int(line.split('-')[0])
					fl_2=int(line.split('-')[1])
					fl_3=int(line.split('-')[2])
					fl_4=int(line.split('-')[3])

					f.write('dihedral: %d,%d,%d,%d,auto\n' % (fl_1,fl_2,fl_3,fl_4))

		f.write('$scan\n')   
		f.write('1: %lf,%lf,%d\n' % (start_angle,start_angle+(n_steps*step_size),n_steps))   

	if gbsa_flag==None:
		os.system('%s %s --opt -c %d -I cntrl_xtb.in > xtb_raw_out' % (xtb_exe,mol_xyz,formal_charge))
	else:
		os.system('%s %s --opt -c %d -I cntrl_xtb.in --alpb %s > xtb_raw_out' % (xtb_exe,mol_xyz,formal_charge,gbsa_flag))

	# get the energies and convert to kcal/mol
	with open('xtbscan.log','r') as f:
		for line in f:
			if 'energy:' in line:
				ene_values.append(float(line.split()[1]))
			elif 'SCF done' in line:
				ene_values.append(float(line.split()[2]))

	zero_point=min(ene_values)
	for i in range(0,len(ene_values)):
		ene_values[i]=(ene_values[i]-zero_point)*627.503

	# convert xyz structures to sdf
	with open('xtbscan.log','r') as f:
		data=f.readlines()
	
	for i in range(0,n_steps):
		local_mol=Chem.Mol(rd_mol)

		local_data=data[(i*len(rd_mol.GetAtoms()))+2+(i*2):(i*len(rd_mol.GetAtoms()))+2+(i*2)+len(rd_mol.GetAtoms())]

		for j in range(0,len(rd_mol.GetAtoms())):
			x=float(local_data[j].split()[1])
			y=float(local_data[j].split()[2])
			z=float(local_data[j].split()[3])

			local_mol.GetConformer().SetAtomPosition(j,(x,y,z))
		
		sd_mols.append(local_mol)
		dihed_values.append(rdMolTransforms.GetDihedralDeg(local_mol.GetConformer(),torsion_list[0]-1,torsion_list[1]-1,torsion_list[2]-1,torsion_list[3]-1))

	output_profile=None
	if output_name is not None:
		output_name=str(output_name)
		output_profile=str(output_name)
	else:
		output_name='torsion_mols'
		output_profile='torsion_profile'

	writer=SDWriter(output_name+'_xtb.sdf')
	counter=0
	for mol in sd_mols:
		this_scan=('%d-%d-%d-%d' % (torsion_list[0],torsion_list[1],torsion_list[2],torsion_list[3]))

		pm=PropertyMol(mol)
		pm.SetProp('xtb-kcal',ene_values[counter])
		pm.SetProp('torsion-angle',dihed_values[counter])
		pm.SetProp('torsion-scanned',this_scan)
		writer.write(pm)
		counter+=1
	writer.flush()

	# output the energy profile
	with open(output_profile+'_xtb.dat','w') as f_out:
		f_out.write('# scan results %d %d %d %d\n' % (torsion_list[0],torsion_list[1],torsion_list[2],torsion_list[3]))
		for i in range(0,len(dihed_values)):
			f_out.write('%lf %lf\n' % (float(dihed_values[i]),float(ene_values[i])))

	os.system('rm %s cntrl_xtb.in charges wbo xtbopt.log xtbscan.log xtbopt.xyz xtbrestart xtb_raw_out .xtboptok xtbtopo.mol' % (mol_xyz))

###############################################################################

if __name__ == '__main__':
	if distutils.spawn.find_executable('xtb'):
		xtb_exe='xtb'
	else:
		print('Please load XTB module\n')
		sys.exit()

	parser = argparse.ArgumentParser(description='XTB dihedral scan.\n')
	parser.add_argument('-i',help='Input mol file (MOL2 or SDF)',required=True)
	parser.add_argument('-n',help='Output name',required=False)
	parser.add_argument('-s',help='Scan dihedral flag',action='store',nargs=4)
	parser.add_argument('-f',help='Dihed fix flag',action='store',nargs=4)
	parser.add_argument('-fl',help='Diheds to fix list',required=False)
	parser.add_argument('-c',help='Formal charge',action='store')
	parser.add_argument('-o',help='Optimize flag',action='store_true')
	parser.add_argument('-gbsa',help='Solvent model for XTB',action='store')

	args=vars(parser.parse_args())

	if mol_type(args['i'])=='sdf':
		suppl=Chem.SDMolSupplier(args['i'],removeHs=False)
		mol=suppl[0]

	elif mol_type(args['i'])=='mol2':
		mol=Chem.MolFromMol2File(args['i'],removeHs=False)

	# Get formal charge from RDKit
	rdchg=int(rdmolops.GetFormalCharge(mol))

	if args['c']!=None:
		formal_charge=int(args['c'])
	else:
		formal_charge=rdchg
	
	if rdchg!=formal_charge:
		print('\n')
		print('Warning: charge specified mis-match',formal_charge,rdchg)
		print('\nExiting\n')
		sys.exit()

	# Make XYZ file
	convert_sdf(mol)

	if (args['o']==False and args['s']!=None):
		# change 36 -> number of torsion scan points
		# change 10 -> degree increment
		print('torsion scan')
		print('\nFormal charge: %d' % (formal_charge))
		xtb_torsion_scan(mol,'tmpAni.xyz',args['s'],args['f'],args['fl'],36,10,xtb_exe,args['gbsa'],formal_charge,args['n'])
	elif (args['o']==True and args['s']==None):
		print('run opt')
		print('\nFormal charge: %d' % (formal_charge))
		xtb_opt(mol,'tmpAni.xyz',xtb_exe,args['gbsa'],formal_charge,args['n'],args['f'],args['fl'])
	else:
		print('Error: args not handled')
		sys.exit(0)

