#!/usr/bin/env python
import sys
import os.path
import math
import numpy as np
import mdtraj as md
import argparse
from scipy.optimize import curve_fit
import time
from example import binning 

################################################################################
# pyWATMD. Run as:
# pyWATMD.py -i input_trajectory -p system_topology -nf frame_process -o output
################################################################################

# Grid spacing
grid_space=1

##### Parse input line
parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="Trajectory for WATMD")
parser.add_argument("-p", type=str, help="Topology file for system")
parser.add_argument("-nf", type=int, help="Number of frames")
parser.add_argument("-o", type=str, help="Output WATMD pdb")

# Required inputs
traj=None
prmtop=None
tot_frames=None
output_pdb=None

# Check command line arguments exist 
if len(sys.argv)==1 or len(sys.argv)==2:
	parser.print_help()
	sys.exit(1)
elif len(sys.argv)>2:
	args = parser.parse_args()
	if (args.i != None and args.nf != None and args.p != None):
		if (os.path.isfile(args.i) and os.path.isfile(args.p)):
			traj=args.i
			prmtop=args.p
			tot_frames=args.nf
		elif (os.path.isfile(args.i) and not os.path.isfile(args.p)):
			print('Cannot find topology file: ',args.p)
			sys.exit(1)
		elif (os.path.isfile(args.p) and not os.path.isfile(args.i)):
			print('Cannot find trajectory file',args.i)
			sys.exit(1)
		elif (not os.path.isfile(args.i) and not os.path.isfile(args.p)):
			print('Cannot find topology or trajectory')
			sys.exit(1)
	else:
		parser.print_help()
		sys.exit(1)

	if args.o != None:
		output_pdb=args.o

		# Print to output
		if len(output_pdb.split('.'))>0:
        		if output_pdb.split('.')[-1]=='pdb':
                		output_pdb=output_pdb.split('.')[0]

	else:
		print('Output WATMD pdb not set: using default watmd_out.pdb \n')
		output_pdb='watmd_out'

if tot_frames!=None and prmtop!=None and traj!=None and output_pdb!=None:
	print('\n')
	print('Processing %s with topology %s and %d frames. Output is %s.pdb\n' % (traj,prmtop,tot_frames,output_pdb))
else:
	print('Error: input options not set')
	sys.exit(1)

##### Define model function to be used to fit to the data
def gauss(x, *p):
	A, mu, sigma = p
	return A*np.exp(-(x-mu)**2/(2.*sigma**2))

################################################################################

# Load toplogy
topology=md.load_prmtop(prmtop)

# Select water atoms 
water_indices=topology.select("water")

# Get min and max coords of water atoms in first frame
with md.open(traj,mode='r') as f:
	for frame in range(0,1):
		xyz,junk_time,cell_lengths,cell_angles=f.read(n_frames=1)
		x_min=np.amin(xyz[0,water_indices[:],0])
		x_max=np.amax(xyz[0,water_indices[:],0])

		y_min=np.amin(xyz[0,water_indices[:],1])
		y_max=np.amax(xyz[0,water_indices[:],1])

		z_min=np.amin(xyz[0,water_indices[:],2])
		z_max=np.amax(xyz[0,water_indices[:],2])

x_box=np.arange(math.floor(x_min),math.ceil(x_max)+grid_space,grid_space)
y_box=np.arange(math.floor(y_min),math.ceil(y_max)+grid_space,grid_space)
z_box=np.arange(math.floor(z_min),math.ceil(z_max)+grid_space,grid_space)

# Print to output
print('Box limits in X: %5.3f -> %5.3f\n' % (math.floor(x_min),math.ceil(x_max)))

#########################################################################

# Print to output
print('Reading frames ... \n')

# Water types: oxy=1, hydr=0
water_type=np.zeros((water_indices.shape[0]))

for i in range(0,water_type.shape[0],3):
	water_type[i]=1

# Init grids
vox_O_count=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))
vox_O_dx=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))
vox_O_dy=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))
vox_O_dz=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))

vox_H_count=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))
vox_H_dx=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))
vox_H_dy=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))
vox_H_dz=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))

start=time.time()

# Read frames
with md.open(traj,mode='r') as f:
	for frame in range(0,tot_frames):
		if frame%1000==0:
			print('Frame: %d' % (frame))

		xyz,junk_time,cell_lengths,cell_angles=f.read(n_frames=1,atom_indices=water_indices[:])

		# Catch error that tot_frames > actual number of frames. Is there a way to catch this earlier ... ?
		try:
			xyz[0,0,0]
		except:
			print('Unable to load frames. nframe > trajectory frame length')
			sys.exit(1)

		count,dx,dy,dz=binning(water_type.astype(np.int32),xyz.astype(np.float64),x_box.astype(np.float64),y_box.astype(np.float64),z_box.astype(np.float64),1)

		vox_O_count=vox_O_count+count
		vox_O_dx=vox_O_dx+dx
		vox_O_dy=vox_O_dy+dy
		vox_O_dz=vox_O_dz+dz

		count,dx,dy,dz=binning(water_type.astype(np.int32),xyz.astype(np.float64),x_box.astype(np.float64),y_box.astype(np.float64),z_box.astype(np.float64),0)

		vox_H_count=vox_H_count+count
		vox_H_dx=vox_H_dx+dx
		vox_H_dy=vox_H_dy+dy
		vox_H_dz=vox_H_dz+dz

end=time.time()
#print('timeit: ',end-start)

# More grids #
vox_bulk=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))

vox_O_meanx=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))
vox_O_meany=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))
vox_O_meanz=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))
vox_O_ocd=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))
vox_O_ohd=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))

vox_H_meanx=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))
vox_H_meany=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))
vox_H_meanz=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))
vox_H_hcd=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))
vox_H_hod=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))

vox_surf=np.zeros((np.shape(x_box)[0],np.shape(y_box)[0],np.shape(z_box)[0]))

vox_O_ocd.fill(-1)
vox_H_hcd.fill(-1)
vox_O_ohd.fill(-1)
vox_H_hod.fill(-1)

# Print to output
print('\n')
print('Fill grid array ... \n')

# Normalise, get ocd, hcd distances 
for i in range(0,np.shape(x_box)[0]):
	for j in range(0,np.shape(y_box)[0]):
		for k in range(0,np.shape(z_box)[0]):
			if vox_O_count[i][j][k]>0:
				vox_O_meanx[i][j][k]=x_box[i] + (vox_O_dx[i][j][k]/vox_O_count[i][j][k])
				vox_O_meany[i][j][k]=y_box[j] + (vox_O_dy[i][j][k]/vox_O_count[i][j][k])
				vox_O_meanz[i][j][k]=z_box[k] + (vox_O_dz[i][j][k]/vox_O_count[i][j][k])
				vox_O_ocd[i][j][k]=math.sqrt((x_box[i]+(float(grid_space)/2) - vox_O_meanx[i][j][k])**2 + (y_box[j]+(float(grid_space)/2) - vox_O_meany[i][j][k])**2 + (z_box[k]+(float(grid_space)/2) - vox_O_meanz[i][j][k])**2 )
				
			if vox_H_count[i][j][k]>0:
				vox_H_meanx[i][j][k]=x_box[i] + (vox_H_dx[i][j][k]/vox_H_count[i][j][k])
				vox_H_meany[i][j][k]=y_box[j] + (vox_H_dy[i][j][k]/vox_H_count[i][j][k])
				vox_H_meanz[i][j][k]=z_box[k] + (vox_H_dz[i][j][k]/vox_H_count[i][j][k])
				vox_H_hcd[i][j][k]=math.sqrt((x_box[i]+(float(grid_space)/2) - vox_H_meanx[i][j][k])**2 + (y_box[j]+(float(grid_space)/2) - vox_H_meany[i][j][k])**2 + (z_box[k]+(float(grid_space)/2) - vox_H_meanz[i][j][k])**2 )

			# ohd, hod distances
			if (vox_O_count[i][j][k]>0 and vox_H_count[i][j][k]>0):
				vox_O_ohd[i][j][k]=math.sqrt((vox_O_meanx[i][j][k]-vox_H_meanx[i][j][k])**2 + (vox_O_meany[i][j][k]-vox_H_meany[i][j][k])**2 + (vox_O_meanz[i][j][k]-vox_H_meanz[i][j][k])**2 )
				vox_H_hod[i][j][k]=math.sqrt((vox_H_meanx[i][j][k]-vox_O_meanx[i][j][k])**2 + (vox_H_meany[i][j][k]-vox_O_meany[i][j][k])**2 + (vox_H_meanz[i][j][k]-vox_O_meanz[i][j][k])**2 )

# ohd, hod distances
for i in range(1,np.shape(x_box)[0]-1,1):
	for j in range(1,np.shape(y_box)[0]-1,1):
		for k in range(1,np.shape(z_box)[0]-1,1):
			if (vox_O_count[i][j][k]>0 and vox_H_count[i][j][k]==0):
				O_ohd_min=999.0
				for i1 in range(i-1,i+2,1):
					for j1 in range(j-1,j+2,1):
						for k1 in range(k-1,k+2,1):
							vox_O_ohd[i][j][k]=math.sqrt((vox_O_meanx[i][j][k]-vox_H_meanx[i1][j1][k1])**2 + (vox_O_meany[i][j][k]-vox_H_meany[i1][j1][k1])**2 + (vox_O_meanz[i][j][k]-vox_H_meanz[i1][j1][k1])**2 )
							O_ohd_min=min(vox_O_ohd[i][j][k],O_ohd_min)
							vox_O_ohd[i][j][k]=O_ohd_min
			elif (vox_O_count[i][j][k]==0 and vox_H_count[i][j][k]>0):
				H_hod_min=999.0
				for i1 in range(i-1,i+2,1):
					for j1 in range(j-1,j+2,1):
						for k1 in range(k-1,k+2,1):
							vox_H_hod[i][j][k]=math.sqrt((vox_H_meanx[i][j][k]-vox_O_meanx[i1][j1][k1])**2 + (vox_H_meany[i][j][k]-vox_O_meany[i1][j1][k1])**2 + (vox_H_meanz[i][j][k]-vox_O_meanz[i1][j1][k1])**2 )
							H_hod_min=min(vox_H_hod[i][j][k],H_hod_min)
							vox_H_hod[i][j][k]=H_hod_min

# Print to output
print('Do gaussian fits ... \n')

# Fitting parameters
bin_width=0.01
guess_mean=0.5
guess_sigma=0.01
guess_coeff=1.0/(math.sqrt(2.0*math.pi*guess_sigma))

##### VOX WAT O OCCUPANCY
hist,bins=np.histogram(vox_O_count,np.arange(0,tot_frames+2,1))
hist[0]=0

guess_O_occ_mean=np.average(vox_O_count,weights=(vox_O_count>0))
#print 'Guess oxygen occupancy: ',guess_O_occ_mean

p0=[1., guess_O_occ_mean, 1.]
try:
        coeff, var_matrix = curve_fit(gauss, bins[:-1], hist, p0=p0)
        vox_O_mean=float(abs(coeff[1]))
        vox_O_std=float(abs(coeff[2]))
        print('vox_O_mean: %lf  vox_O_std: %lf' % (vox_O_mean,vox_O_std))
except:
        print('vox oxygen occupancy gauss fit failed')
        sys.exit(1)
        pass

##### VOX WAT H OCCUPANCY
hist,bins=np.histogram(vox_H_count,np.arange(0,2*tot_frames+2,1))
hist[0]=0

guess_H_occ_mean=np.average(vox_H_count,weights=(vox_H_count>0))
#print 'Guess hydrogen occupancy: ',guess_H_occ_mean

p0=[1., guess_H_occ_mean, 1.]
try:
        coeff, var_matrix = curve_fit(gauss, bins[:-1], hist, p0=p0)
        vox_H_mean=abs(coeff[1])
        vox_H_std=abs(coeff[2])
        print('vox_H_mean: %lf  vox_H_std: %lf' % (vox_H_mean,vox_H_std))
except:
        print('vox hydrogen occupancy gauss fit failed')
        sys.exit(1)
        pass

##### VOX WAT O CENTROID
hist,bins=np.histogram(vox_O_ocd,np.arange(0,1,bin_width))

p0=[guess_coeff, guess_mean, guess_sigma]
try:
        coeff, var_matrix = curve_fit(gauss, bins[:-1], hist, p0=p0)
        vox_ocd_mean=abs(coeff[1])
        vox_ocd_std=abs(coeff[2])
        print('vox_ocd_mean: %lf  vox_ocd_std: %lf' % (vox_ocd_mean,vox_ocd_std))
except:
        print('vox oxygen centroid gauss fit failed')
        sys.exit(1)
        pass

##### VOX WAT H CENTROID
hist,bins=np.histogram(vox_H_hcd,np.arange(0,1,bin_width))

p0=[guess_coeff, guess_mean, guess_sigma]
try:
        coeff, var_matrix = curve_fit(gauss, bins[:-1], hist, p0=p0)
        vox_hcd_mean=abs(coeff[1])
        vox_hcd_std=abs(coeff[2])
        print('vox_hcd_mean: %lf  vox_hcd_std: %lf' % (vox_hcd_mean,vox_hcd_std))
except:
        print('vox hydrogen centroid gauss fit failed')
        sys.exit(1)
        pass

##### VOX WAT O -> H
hist,bins=np.histogram(vox_O_ohd,np.arange(0,2,bin_width))

p0=[guess_coeff, guess_mean, guess_sigma]
try:
        coeff, var_matrix = curve_fit(gauss, bins[:-1], hist, p0=p0)
        vox_ohd_mean=abs(coeff[1])
        vox_ohd_std=abs(coeff[2])
        print('vox_ohd_mean: %lf  vox_ohd_std: %lf' % (vox_ohd_mean,vox_ohd_std))
except:
        print('vox oxygen -> H gauss fit failed')
        sys.exit(1)
        pass

##### VOX WAT H -> O
hist,bins=np.histogram(vox_H_hod,np.arange(0,2,bin_width))

p0=[guess_coeff, guess_mean, guess_sigma]
try:
        coeff, var_matrix = curve_fit(gauss, bins[:-1], hist, p0=p0)
        vox_hod_mean=abs(coeff[1])
        vox_hod_std=abs(coeff[2])
        print('vox_hod_mean: %lf  vox_hod_std: %lf' % (vox_hod_mean,vox_hod_std))
except:
        print('vox hydrogen -> O gauss fit failed')
        sys.exit(1)
        pass

# surface - 5A buffer from box edge
print('\nSurface calculation ...\n')
buffer_box=5.0*(1.0/float(grid_space))

for i in range(1,np.shape(x_box)[0]-1,1):
	for j in range(1,np.shape(y_box)[0]-1,1):
		for k in range(1,np.shape(z_box)[0]-1,1):
			if (vox_O_count[i][j][k]>(1.0/3.0)*vox_O_mean) or (vox_H_count[i][j][k]>(1.0/3.0)*vox_H_mean):
				try:
					x_fill_min=np.nonzero(vox_O_count[:,j,k])[0][0]
					x_fill_max=np.nonzero(vox_O_count[:,j,k])[0][-1]
					y_fill_min=np.nonzero(vox_O_count[i,:,k])[0][0]
					y_fill_max=np.nonzero(vox_O_count[i,:,k])[0][-1]
					z_fill_min=np.nonzero(vox_O_count[i,j,:])[0][0]
					z_fill_max=np.nonzero(vox_O_count[i,j,:])[0][-1]
				except:
					x_fill_min=100
					x_fill_max=-1
					y_fill_min=100
					y_fill_max=-1
					z_fill_min=100
					z_fill_max=-1

				for i1 in range(i-1,i+2,1):
					for j1 in range(j-1,j+2,1):
						for k1 in range(k-1,k+2,1):
							if (vox_O_count[i1][j1][k1]<(1.0/3.0)*vox_O_mean) or (vox_H_count[i1][j1][k1]<(1.0/3.0)*vox_H_mean):

								if (x_fill_min+buffer_box<i<x_fill_max-buffer_box) and (y_fill_min+buffer_box<j<y_fill_max-buffer_box) and (z_fill_min+buffer_box<k<z_fill_max-buffer_box):
									vox_surf[i][j][k]=1

# Print to output
print('Finding bulk water ... \n')

n_sigma_bulk=3
for i in range(0,np.shape(x_box)[0]):
	for j in range(0,np.shape(y_box)[0]):
		for k in range(0,np.shape(z_box)[0]):
			# Check cells for bulk water
			bulk_condition=0
			if (vox_O_ocd[i][j][k]>(vox_ocd_mean-n_sigma_bulk*vox_ocd_std) and vox_O_ocd[i][j][k]<(vox_ocd_mean+n_sigma_bulk*vox_ocd_std)):
				bulk_condition+=1
			if (vox_H_hcd[i][j][k]>(vox_hcd_mean-n_sigma_bulk*vox_hcd_std) and vox_H_hcd[i][j][k]<(vox_hcd_mean+n_sigma_bulk*vox_hcd_std)):
				bulk_condition+=1
			if (vox_O_ohd[i][j][k]>(vox_ohd_mean-n_sigma_bulk*vox_ohd_std) and vox_O_ohd[i][j][k]<(vox_ohd_mean+n_sigma_bulk*vox_ohd_std)):
				bulk_condition+=1
			if (vox_H_hod[i][j][k]>(vox_hod_mean-n_sigma_bulk*vox_hod_std) and vox_H_hod[i][j][k]<(vox_hod_mean+n_sigma_bulk*vox_hod_std)):
				bulk_condition+=1
			if (vox_O_count[i][j][k]>(vox_O_mean-n_sigma_bulk*vox_O_std) and vox_O_count[i][j][k]<(vox_O_mean+n_sigma_bulk*vox_O_std)) and vox_O_count[i][j][k]>0:
				bulk_condition+=1
			if bulk_condition==6:
				vox_bulk[i][j][k]=-1

# Print to output
print('Writing output: %s.pdb\n' % (output_pdb))

# Write output pdb results
fout=open(output_pdb+'.pdb','w')
fout.write('HEADER %s_CVX\n' % (output_pdb))
atm_count=1
for i in range(0,np.shape(x_box)[0]):
	for j in range(0,np.shape(y_box)[0]):
		for k in range(0,np.shape(z_box)[0]):	
			if vox_bulk[i][j][k]>-1:
				if (vox_O_count[i][j][k]>0 or vox_H_count[i][j][k]>0):
					if vox_O_count[i][j][k]>0:
						vox_O_meanx[i][j][k]=(vox_O_meanx[i][j][k]-x_box[i])
						vox_O_meany[i][j][k]=(vox_O_meany[i][j][k]-y_box[j])
						vox_O_meanz[i][j][k]=(vox_O_meanz[i][j][k]-z_box[k])

					if vox_H_count[i][j][k]>0:
						vox_H_meanx[i][j][k]=(vox_H_meanx[i][j][k]-x_box[i])
						vox_H_meany[i][j][k]=(vox_H_meany[i][j][k]-y_box[j])
						vox_H_meanz[i][j][k]=(vox_H_meanz[i][j][k]-z_box[k])

					if (vox_O_count[i][j][k]>0 and vox_H_count[i][j][k]>0):
						VoxChg = ((0.5*vox_H_count[i][j][k]) - vox_O_count[i][j][k])/(vox_O_count[i][j][k]+vox_H_count[i][j][k])
						NumWat = ((vox_O_count[i][j][k]+(0.5*vox_H_count[i][j][k]))/(vox_O_mean+(0.5*vox_H_mean)))
						dx=((x_box[i]+vox_H_meanx[i][j][k]) - (x_box[i]+vox_O_meanx[i][j][k]))
						Xo=(x_box[i] + vox_O_meanx[i][j][k]) + dx*(vox_H_count[i][j][k]/(vox_O_count[i][j][k]+vox_H_count[i][j][k]))
						dy=((y_box[j]+vox_H_meany[i][j][k]) - (y_box[j]+vox_O_meany[i][j][k]))
						Yo=(y_box[j] + vox_O_meany[i][j][k]) + dy*(vox_H_count[i][j][k]/(vox_O_count[i][j][k]+vox_H_count[i][j][k]))
						dz=((z_box[k]+vox_H_meanz[i][j][k]) - (z_box[k]+vox_O_meanz[i][j][k]))
						Zo=(z_box[k] + vox_O_meanz[i][j][k]) + dz*(vox_H_count[i][j][k]/(vox_O_count[i][j][k]+vox_H_count[i][j][k]))
					elif (vox_O_count[i][j][k]>0 and vox_H_count[i][j][k]==0):
						VoxChg=-1
						NumWat=(vox_O_count[i][j][k]/vox_O_mean)
						Xo=(x_box[i] + vox_O_meanx[i][j][k])
						Yo=(y_box[j] + vox_O_meany[i][j][k])
						Zo=(z_box[k] + vox_O_meanz[i][j][k])
					elif (vox_O_count[i][j][k]==0 and vox_H_count[i][j][k]>0):
						VoxChg=0.5
						NumWat=(vox_H_count[i][j][k]/vox_H_mean)
						Xo=(x_box[i] + vox_H_meanx[i][j][k])
						Yo=(y_box[j] + vox_H_meany[i][j][k])
						Zo=(z_box[k] + vox_H_meanz[i][j][k])
					# Hydrogen charge double 
					if VoxChg>0:
						VoxChg=VoxChg*2
					# skip printing of atm_count due to overflow
					line=str("ATOM%7d   Ne CVX     1    %8.3f%8.3f%8.3f %5.2f %5.2f\n" % (int(1),Xo,Yo,Zo,VoxChg,NumWat))
					atm_count+=1
					fout.write(line)

fout.write('HEADER %s_SUR\n' % (output_pdb))
atm_count=1
for i in range(0,np.shape(x_box)[0]):
	for j in range(0,np.shape(y_box)[0]):
		for k in range(0,np.shape(z_box)[0]):
			if (vox_surf[i][j][k])>0:
				NumWat=(vox_O_count[i][j][k]/vox_O_mean)
				Xo=(x_box[i] + vox_O_meanx[i][j][k])
				Yo=(y_box[j] + vox_O_meany[i][j][k])
				Zo=(z_box[k] + vox_O_meanz[i][j][k])
				line=str("ATOM%7d   He SUR     2    %8.3f%8.3f%8.3f %5.2f %5.2f\n" % (int(1),Xo,Yo,Zo,float(1.0),NumWat))
				atm_count+=1
				fout.write(line)

fout.write('END\n')

fout.close()

