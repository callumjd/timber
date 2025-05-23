#!/usr/bin/env python
import math
import numpy as np
import argparse

#
# SPLIT method to determine counter-ions
# Machado, Pantano JCTC 2020, 16, 3
# https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953
#

def return_salt(nwat,conc,charge):

	N0=(nwat*conc)/56.0

	Npos=int(math.ceil(N0-(charge/2)))
	Nneg=int(math.ceil(N0+(charge/2)))

	return Npos,Nneg


if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Salt concentration solver\n')

	parser.add_argument('-nw',help='Number of waters',type=int,required=True)
	parser.add_argument('-conc',help='Salt concentration (M)',type=float,required=False,default=0.15)
	parser.add_argument('-q',help='System charge',type=int,required=False,default=0)

	args=vars(parser.parse_args())

	nwat=int(args['nw'])
	conc=float(args['conc'])
	charge=float(args['q'])

	Npos,Nneg=return_salt(nwat,conc,charge)

	print('\naddionsrand mol Na+ %d Cl- %d\n' % (Npos,Nneg))

