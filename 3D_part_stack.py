# This program takes 3D particle data of 100 halo sample from Gerard Lemson from the MDB
# and stacks the data by mass bin and uses the M_Phi technique. Note:(each bin is an ensemble cluster)

# last update: 1/29/13

##########

from numpy import *
from stack_class_3D import *
from astStats import *
from flux_caustics_ideal import *
from matplotlib.pyplot import *
import cosmolopy.distance as cd

## DEFINE CONSTANTS ##

h = 0.72 		# Hubble Constant / 100.0
r_limit = 2		# Radius Limit of data in Mpc
H0 = h*100.0		# Hubble constant
q = 10.0
c = 300000.0
cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'h':H0/100.0}
cosmo = cd.set_omega_k_0(cosmo)
halo_num = 100		# Total number of halos
run_num = [0,10]	# Number of halos to run program over, must be continuous
bin_range = 10		# Number of halos per bin
gal_num = 1000		# Number of particles stacked per halo for ensemble clusters

## DEFINE FLAGS ##

use_mems = False
use_vdisp = True

## INITIALIZATION ##

G = galaxies()
P = particles()
C = caustic()
U = universal()

### PROGRAM ###

print '...loading halos'

HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.load_halos(h)

HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.sort_halos(HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z)

print '...loading particles'

R, V, PPX, PPY, PPZ = P.configure_particles(HaloID,h,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_limit,R_crit200,HVD,halo_num,gal_num,run_num)

print '...binning data'
# All variables beginning with 'ENC_' stand for ensemble cluster, same for *bin variables

ENC_R,ENC_V,ENC_M200,ENC_R200,ENC_HVD,ENC_SRAD,ENC_ESRAD = P.bin_data(HaloID,R,V,SRAD,ESRAD,M_crit200,R_crit200,HVD,halo_num,bin_range,run_num,gal_num)

print '...running caustic'

x_range,ENC_INF_NFWMASS,ENC_DIA_NFWMASS,ENC_INF_CAUMASS,ENC_DIA_CAUMASS,ENC_INF_MPROF,ENC_INF_NFW,ENC_INF_CAU,ENC_DIA_MPROF,ENC_DIA_NFW,ENC_DIA_CAU = P.kernel_caustic_masscalc(ENC_R,ENC_V,ENC_M200,ENC_R200,ENC_SRAD,ENC_ESRAD,ENC_HVD,halo_num,bin_range,gal_num,H0,q,r_limit,run_num,use_mems)



