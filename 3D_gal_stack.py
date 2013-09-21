# This program takes 3D galaxy data of 100 halo sample from Gerard Lemson from the MDB
# and stacks the data by mass bin and uses the M_Phi technique. Note:(each bin is an ensemble cluster)

# last update: 1/24/13

######

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
line_num = 1			# Number of iterations (lines of sight, however b/c 3D program, only controls iterations)
halo_num = 100			# Total number of halos
gal_num = 10			# Number of galaxies stacked per halo for en. clusters
bin_range = 10			# Number of halos per ensemble cluster
run_num = [0,10]		# Number of ensemble clusters to run caustic mass calc over

## DEFINE FLAGS ##

use_mems = False
use_vdisp = True

## INITIALIZATION ##

G = galaxies()
C = caustic()
U = universal()

### PROGRAM ###

print '...loading halos'

HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.load_halos(h)

HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.sort_halos(HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z)

print '...loading galaxies'

R, V, MAGS, GPX, GPY, GPZ = G.configure_galaxies(HaloID,h,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_limit,R_crit200,HVD,halo_num)

print '...binning data'	
# All variables beginning with 'ENC_' stand for ensemble cluster variable, similar to *bin variables

ENC_R,ENC_V,MAGbin,ENC_M200,ENC_R200,ENC_HVD,ENC_SRAD,ENC_ESRAD,ENC_GPX,ENC_GPY,ENC_GPZ = G.bin_data_mag(HaloID,R,V,MAGS,SRAD,ESRAD,M_crit200,R_crit200,HVD,halo_num,bin_range,gal_num,GPX,GPY,GPZ)

print '...running caustic'

x_range,ENC_INF_NFWMASS,ENC_DIA_NFWMASS,ENC_INF_CAUMASS,ENC_DIA_CAUMASS,ENC_INF_MPROF,ENC_INF_NFW,ENC_INF_CAU,ENC_DIA_MPROF,ENC_DIA_NFW,ENC_DIA_CAU = G.kernel_caustic_masscalc(ENC_R,ENC_V,ENC_M200,ENC_R200,ENC_SRAD,ENC_ESRAD,ENC_HVD,halo_num,bin_range,gal_num,H0,q,r_limit,run_num,use_mems)

print ''
bias1 = mean( (ENC_M200-ENC_INF_NFWMASS) / ENC_M200 )
bias2 = mean( abs(ENC_M200-ENC_INF_NFWMASS) / ENC_M200 )
bias3 = mean( log(ENC_INF_NFWMASS/ENC_M200) )


