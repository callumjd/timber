cimport cython
import numpy as np
cimport numpy as np

def find_bin(double x, double[:] r):
	cdef int N = r.shape[0]
	cdef int i
	cdef int b

	for i in range(N):
		if x==r[i]:
			b=i
		elif r[i]<x<r[i+1]:
			b=i
		
	return b

def binning(int[:] N, double[:, :, :] xyz, double[:] X, double[:] Y, double[:] Z, int type): 
	cdef double[:, :, :] vox_count=np.zeros((X.shape[0],Y.shape[0],Z.shape[0]))
	cdef double[:, :, :] vox_dx=np.zeros((X.shape[0],Y.shape[0],Z.shape[0]))
	cdef double[:, :, :] vox_dy=np.zeros((X.shape[0],Y.shape[0],Z.shape[0]))
	cdef double[:, :, :] vox_dz=np.zeros((X.shape[0],Y.shape[0],Z.shape[0]))

	cdef double x_wat, y_wat, z_wat
	cdef int i
	cdef int xb,yb,zb
	cdef double dx, dy, dz

	for i in range(N.shape[0]):
		if N[i]==type:
			x_wat=xyz[0,i,0]
			y_wat=xyz[0,i,1]
			z_wat=xyz[0,i,2]

			if (x_wat>X[0] and x_wat<X[-1]):
				xb=find_bin(x_wat,X)
				dx=abs(x_wat-X[xb])

			if (y_wat>Y[0] and y_wat<Y[-1]):
				yb=find_bin(y_wat,Y)
				dy=abs(y_wat-Y[yb])

			if (z_wat>Z[0] and z_wat<Z[-1]):
				zb=find_bin(z_wat,Z)
				dz=abs(z_wat-Z[zb])

			if ((-1<xb<X.shape[0]+1) and (-1<yb<Y.shape[0]+1) and (-1<zb<Z.shape[0]+1)):
				vox_count[xb][yb][zb]+=1
				vox_dx[xb][yb][zb]+=dx
				vox_dy[xb][yb][zb]+=dy
				vox_dz[xb][yb][zb]+=dz

	return vox_count,vox_dx,vox_dy,vox_dz	
