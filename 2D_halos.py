## This program uses the Caustic and M_PHI technique to produce mass estimates for Millennium Database Halos

###########

from numpy import *
from stack_class_2D import *
from astStats import biweightScale
from flux_caustics_nideal import *
from matplotlib.pyplot import *
from numpy.random import randint

## DEFINE FLAGS ##

use_mems = False
use_vdisp = True	# Feed table value hvd
use_gals = True		# Use galaxies or particles?

## DEFINE CONSTANTS ##

h = 0.72 			# Hubble Constant / 100.0
r_limit = 2			# Radius Limit of data in R_200
H0 = h*100.0			# Hubble constant
q = 10.0
c = 300000.0
bin_range = 1				# Needed b/c technically it is working on ensemble code
halo_num = 100				# Number of halos in sample
run_num = [0,100]			# Number of halos to run program over, particles
vlimit = 3500.0	
gal_num = 100

	
## INITIALIZATION ##

C = caustic()
G = galaxies()
U = universal()

## DEFINE FUNCTIONS ##

def get_galsbigGuo(self,ID,H0,Z):
	fileID = np.loadtxt('/n/Christoq1/MILLENNIUM/particles/cmiller.csv',dtype='string',delimiter=',',skiprows=1,usecols=(0,),unpack=True)
	Nid = np.where(ID==fileID)[0][0]
	fits = pyfits.open('/n/Christoq1/MILLENNIUM/particles/t_'+str(Nid)+'_cmiller_guo.fits')
	xpos,ypos,zpos = fits[1].data.field('x')/(1+Z),fits[1].data.field('y')/(1+Z),fits[1].data.field('z')/(1+Z)
	return (fits[1].data.field('HALOID'),fits[1].data.field('uDUST'),fits[1].data.field('gDUST'),fits[1].data.field('rDUST'),fits[1].data.field('iDUST'),fits[1].data.field('zDUST'),xpos/(H0/100.0),ypos/(H0/100.0),zpos/(H0/100.0),fits[1].data.field('velX'),fits[1].data.field('velY'),fits[1].data.field('velZ'),fits[1].data.field('vvir'))


def get_galsbigParticles(self,ID,H0,Z):
	fileID = np.loadtxt('/n/Christoq1/MILLENNIUM/particles/cmiller.csv',dtype='string',delimiter=',',skiprows=1,usecols=(0,),unpack=True)
	snapnum = np.loadtxt('/n/Christoq1/MILLENNIUM/particles/cmiller.csv',dtype='float',delimiter=',',skiprows=1,usecols=(1,),unpack=True)
	Nid = np.where(ID==fileID)[0][0]
	fits = pyfits.open('/n/Christoq1/MILLENNIUM/particles/t'+str(Nid)+'_cmiller.dat.fits')
	part_x = fits[1].data.field('PPX')/(1+Z)
	part_y = fits[1].data.field('PPY')/(1+Z)
	part_z = fits[1].data.field('PPZ')/(1+Z)
	part_vx = fits[1].data.field('VVX')/np.sqrt(1+Z)
	part_vy = fits[1].data.field('VVY')/np.sqrt(1+Z)
	part_vz = fits[1].data.field('VVZ')/np.sqrt(1+Z)
	return (fits[1].data.field('ID'),np.random.random(part_x.size),np.random.random(part_x.size),np.random.random(part_x.size),np.random.random(part_x.size),np.random.random(part_x.size),part_x/(H0/100.0),part_y/(H0/100.0),part_z/(H0/100.0),part_vx,part_vy,part_vz)



#####	PROGRAM    ######

print '...loading halos'
HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.load_halos(h)
HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.sort_halos(HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z)


print '...loading galaxies'
if use_gals==True:
	R, V, MAGS, GAL_VDISP, GPX, GPY = G.configure_galaxies(HaloID,h,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_limit,vlimit,R_crit200,HVD,halo_num,gal_num)


print '...caustic!'
if use_gals==True:
	x_range,INF_NFWMASS,DIA_NFWMASS,INF_CAUMASS,DIA_CAUMASS,INF_MPROF,INF_NFW,INF_CAU,DIA_MPROF,DIA_NFW,DIA_CAU = G.kernel_caustic_masscalc(R,V,M_crit200,R_crit200,SRAD,ESRAD,HVD,GAL_VDISP,halo_num,bin_range,gal_num,H0,q,r_limit,run_num,use_vdisp,use_mems)












