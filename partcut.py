#This code takes particle data from /n/Christoq1/MILLENNIUM/particles, and cuts the particle data by a general cut, so that when I upload them, I'm not dealing with millions of them.

import pyfits
from numpy import *
from stack_class_3D import *


h = 0.72

U=universal()

HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.load_halos(h)  

HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.sort_halos(HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z)

k=0
for haloid in HaloID:
	id = loadtxt('/n/Christoq1/MILLENNIUM/particles/cmiller.csv', dtype='str', delimiter=',', usecols=(0,), unpack=True)
	id = delete(id,0)	
	index = where(id==str(haloid))
	p = pyfits.open('/n/Christoq1/MILLENNIUM/particles/t'+str(index[0][0])+'_cmiller.dat.fits')	









