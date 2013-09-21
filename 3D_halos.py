# 3d galaxy and particle m_phi code for non-stacked halos.

# last update: 1/24/13

#######

from numpy import *
from stack_class_3D import *
from astStats import *
from flux_caustics_ideal import *
from matplotlib.pyplot import *
import cosmolopy.distance as cd
from numpy.random import randint

## DEFINE FLAGS ##

use_mems = False
use_vdisp = True

use_gals = True			# Use galaxies, or else particles 

## DEFINE CONSTANTS ##

h = 0.72 			# Hubble Constant / 100.0
r_limit = 2			# Radius Limit of data in Mpc
H0 = h*100.0			# Hubble constant
q = 10.0
c = 300000.0
cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'h':H0/100.0}
cosmo = cd.set_omega_k_0(cosmo)
bin_range = 1			# Needed b/c technically it is working on ensemble code
halo_num = 100			# Number of halos in sample
run_num = [0,1]			# Number of halos to run program over, particles

if use_gals==True:
	gal_num = 10		# Number of galaxies per mass calculation
else:
	gal_num = 10000		# Number of particles per mass calculation

## DEFINE FUNCTIONS ##

def load_gals(h,gal_num,HaloID,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_crit200):
	R = []
	V = []
	MAGS = []
	GPX = []
	GPY = []
	GPZ = []	

	ID = loadtxt('/n/Christoq1/nkern/Documents/MDB_milliMil_halodata/Caustic/cmiller.csv',delimiter=',',dtype='str',usecols=(0,),unpack=True)		
	for haloid,k in zip(list(HaloID),list(range(halo_num))):
		IDmatch = where(ID==haloid)[0][0]
		f = pyfits.open('/n/Christoq1/MILLENNIUM/particles/t_'+str(IDmatch)+'_cmiller_guo.fits')
		data = f[1].data
		z,gpx,gpy,gpz,gvx,gvy,gvz,mags = data.field(13),data.field(17),data.field(18),data.field(19),data.field(20),data.field(21),data.field(22),data.field(63)

		gpx,gpy,gpz = ( gpx/(1+Z[k])/h - HPX[k] ),( gpy/(1+Z[k])/h - HPY[k] ),( gpz/(1+Z[k])/h - HPZ[k] )
		gvx,gvy,gvz = gvx-HVX[k],gvy-HVY[k],gvz-HVZ[k]
	
		r = sqrt( (gpx)**2 + (gpy)**2 + (gpz)**2 ) 
		v = sqrt( (gvx)**2 + (gvy)**2 + (gvz)**2 )

		sort = argsort(mags)		# mag sorted
		mags = array(mags[sort])
		r = array(r[sort])		 
		v = array(v[sort])
		gpx,gpy,gpz,gvx,gvy,gvz = array(gpx[sort]),array(gpy[sort]),array(gpz[sort]),array(gvx[sort]),array(gvy[sort]),array(gvz[sort])
		## LIMIT DATA ##
		cut = where((r<=r_limit*R_crit200[k]) & (v<=5000.0) & (v!=0))[0][0:gal_num]
		r,v,mags = r[cut],v[cut],mags[cut]
		gpx,gpy,gpz,gvx,gvy,gvz = gpx[cut],gpy[cut],gpz[cut],gvx[cut],gvy[cut],gvz[cut]	

		R.append(r)
		V.append(v)
		MAGS.append(mags)
		GPX.append(gpx)
		GPY.append(gpy)
		GPZ.append(gpz)
	R = array(R)
	V = array(V)
	MAGS = array(MAGS)
	GPX = array(GPX)
	GPY = array(GPY)
	GPZ = array(GPZ)
	return R,V,MAGS,GPX,GPY,GPZ	

def load_parts(h,HaloID,HVX,HVY,HVZ,Z,run_num):
	R = []
	V = []
	PPX = []
	PPY = []
	PPZ = []
	for haloid,k in zip(list(HaloID[run_num[0]:run_num[1]]),list(arange(run_num[0],run_num[1]))):
		id = loadtxt('/n/Christoq1/MILLENNIUM/particles/cmiller.csv', dtype='str', delimiter=',', usecols=(0,), unpack=True)
		id = delete(id,0)		
		index = where(id==str(haloid))
		p = pyfits.open('/n/Christoq1/MILLENNIUM/particles/t'+str(index[0][0])+'_cmiller.dat.fits')
		data = p[1].data
		ppx = data.field(1)/h/(1+Z[k])
		ppy = data.field(2)/h/(1+Z[k])
		ppz = data.field(3)/h/(1+Z[k])
		pvx = data.field(4)/sqrt(1+Z[k])
		pvy = data.field(5)/sqrt(1+Z[k])
		pvz = data.field(6)/sqrt(1+Z[k])
	
		pvx,pvy,pvy = pvx-HVX[k],pvy-HVY[k],pvz-HVZ[k]
	
		r = sqrt( (ppx**2) + (ppy**2) + (ppz**2) ) 
		v = sqrt( (pvx)**2 + (pvy)**2 + (pvz)**2 )

		r = array(r)
		v = array(v)

		## LIMIT DATA ##
		cut = where((r<=r_limit) & (v<=5000.0))
		r = r[cut]
		v = v[cut]
		ppx,ppy,ppz = ppx[cut],ppy[cut],ppz[cut]
		pick = randint(0,r.size,gal_num)
		r = r[pick]
		v = v[pick]
		ppx,ppy,ppz = ppx[pick],ppy[pick],ppz[pick]
		R.append(r)
		V.append(v)
		PPX.append(ppx)
		PPY.append(ppy)
		PPZ.append(ppz)
		print 'done loading halo',k

	R = array(R)
	V = array(V)
	PPX = array(PPX)
	PPY = array(PPY)
	PPZ = array(PPZ)
	return R, V, PPX, PPY, PPZ
		
## INITIALIZATION ##

U = universal()
P = particles()
G = galaxies()
C = caustic()

### PROGRAM ###

print '...loading halos'

HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.load_halos(h)

HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.sort_halos(HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z)

if use_gals == True:
	print '...loading gals'
	R, V, MAGS, GPX, GPY, GPZ = load_gals(h,gal_num,HaloID,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,R_crit200)

else:
	print '...loading particles'
	R, V, PPX, PPY = load_parts(h,HaloID,HVX,HVY,HVZ,Z,run_num)

print '...caustic!'
if use_gals == True:
	x_range,ENC_INF_NFWMASS,ENC_DIA_NFWMASS,ENC_INF_CAUMASS,ENC_DIA_CAUMASS,ENC_INF_MPROF,ENC_INF_NFW,ENC_INF_CAU,ENC_DIA_MPROF,ENC_DIA_NFW,ENC_DIA_CAU = G.kernel_caustic_masscalc(R,V,M_crit200,R_crit200,SRAD,ESRAD,HVD,halo_num,bin_range,gal_num,H0,q,r_limit,run_num,use_mems)

else:
	x_range,ENC_INF_NFWMASS,ENC_DIA_NFWMASS,ENC_INF_CAUMASS,ENC_DIA_CAUMASS,ENC_INF_MPROF,ENC_INF_NFW,ENC_INF_CAU,ENC_DIA_MPROF,ENC_DIA_NFW,ENC_DIA_CAU = P.kernel_caustic_masscalc(R,V,M_crit200,R_crit200,SRAD,ESRAD,HVD,halo_num,bin_range,gal_num,H0,q,r_limit,run_num,use_mems)	

