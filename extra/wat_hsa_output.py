#!/usr/bin/env python
import argparse
import numpy as np

E_bulk=-9.53 # tip3p bulk water free energy

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='simple parser of HSA output\n')
	parser.add_argument('-i',help='HSA summary txt file',required=True)
	args=vars(parser.parse_args())

	hsa_summary=args['i']

	# Load data, skip solute_acceptors solute_donors string fields
	data=np.genfromtxt(hsa_summary,skip_header=1,usecols=range(0,27))

	with open('clustercenterfile.pdb','r') as f:
		f.readline()
		pdb_clus=f.readlines()

	print('Writing E_tot cluster file: E_rel-clustercenterfile.pdb .. \n')
	# E_tot corrected for bulk
	with open('E_rel-clustercenterfile.pdb','w') as f_out:
		for i in range(0,np.shape(data)[0]):
			str_out=pdb_clus[i][0:54]
			f_out.write(('%s %5.2f %5.2f\n') % (str_out,float(data[i][12]-E_bulk),float(data[i][5])))

	print('Writing G_tot cluster file: G_rel-clustercenterfile.pdb .. \n')
	# Free energy corrected for bulk 
	with open('G_rel-clustercenterfile.pdb','w') as f_out:
		for i in range(0,np.shape(data)[0]):
			str_out=pdb_clus[i][0:54]
			f_out.write(('%s %5.2f %5.2f\n') % (str_out,float((data[i][12]-E_bulk)-data[i][16]),float(data[i][5])))

	print('Writing NBR cluster file: NBR-clustercenterfile.pdb .. \n')
	# Ewwnbr 
	with open('NBR-clustercenterfile.pdb','w') as f_out:
		for i in range(0,np.shape(data)[0]):
			str_out=pdb_clus[i][0:54]
			f_out.write(('%s %5.2f %5.2f\n') % (str_out,float(data[i][13]),float(data[i][5])))

	print('Writing Nhbtot cluster file: HBtot-clustercenterfile.pdb .. \n')
	# Nhbtot
	with open('HBtot-clustercenterfile.pdb','w') as f_out:
		for i in range(0,np.shape(data)[0]):
			str_out=pdb_clus[i][0:54]
			f_out.write(('%s %5.2f %5.2f\n') % (str_out,float(data[i][20]),float(data[i][5])))

