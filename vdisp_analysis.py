
######
from numpy import *
from stack_class_2D import *
from flux_caustics_nideal import *
from matplotlib.pyplot import *
import sys
import random
 
## DEFINE CONSTANTS ##
h = 0.72 		# Hubble Constant / 100.0
r_limit = 1.25		# Radius Limit of data in R_crit200 Scale Factor
vlimit = 3500.0		# Velocity Disp limit in km/s
H0 = h*100.0		# Hubble constant
q = 10.0		# Normalization...
c = 300000.0		# Speed of Light
halo_num = 100		# Total number of halos loaded
beta = 0.2				# Beta constant
gal_num = int(sys.argv[2])		# Number of galaxies stacked per halo for en. clusters
line_num = int(sys.argv[3])		# Number of lines of sight to stack over

run_loc = 'ss_run'+str(sys.argv[4])+''	# What directory do you want to write files in?

ens_num = [int(sys.argv[1])]	#Ensemble Number, for pbs submission w/ job array

## DEFINE FLAGS ##
use_mems = False	# Membership not known a priori
use_vdisp = True	# If false, do clipping
scale_data = False	# Scale data before stacking

use_flux = True		# Using flux or sophie?
write_data = False	# Write out data?

if use_flux == True:
	root = str('/nfs/christoq_ls')
else:
	root = str('/n/Christoq1')

## INITIALIZATION ##
G = galaxies()
C = caustic()
U = universal()

### PROGRAM ###

print '...loading halos'
HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.load_halos(h,root)
HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.sort_halos(HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z)

print '...loading galaxies'
Halo_P,Halo_V,Gal_P,Gal_V,MAGS = G.configure_galaxies(HaloID,h,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_limit,vlimit,R_crit200,HVD,halo_num,gal_num,root)

# Data Arrays, Method 0, Method 1, Method 2
# Method 0: Top Brightest
# Method 1: Random Top Brightest
# Method 2: Ordered sets of Top Brightest

ENC_R0,ENC_V0,ENC_MAG0 = [],[],[]
ENC_R1,ENC_V1,ENC_MAG1 = [],[],[]
ENC_R2,ENC_V2,ENC_MAG2 = [],[],[]
LINE_R0,LINE_V0,LINE_MAG0 = [],[],[]
LINE_R1,LINE_V1,LINE_MAG1 = [],[],[]
LINE_R2,LINE_V2,LINE_MAG2 = [],[],[]
LENGTH = []
##


# Loop over Halo
for k in ens_num:
	enc_r0,enc_v0,enc_mag0 = [],[],[]
	enc_r1,enc_v1,enc_mag1 = [],[],[]
	enc_r2,enc_v2,enc_mag2 = [],[],[]
	line_r0,line_v0,line_mag0 = [],[],[]
	line_r1,line_v1,line_mag1 = [],[],[]
	line_r2,line_v2,line_mag2 = [],[],[]
	length = []

	# Loop over Lines of Sight
	for l in range(line_num):
		# Do Line of Sight
		r,v = U.line_of_sight(Gal_P[k],Gal_V[k],Halo_P[k],Halo_V[k],H0,c)
		mags = MAGS[k]	
		# Sort Gals
		sort = argsort(mags)
		r,v,mags = r[sort],v[sort],mags[sort]
		# Limit Gals
		limit = where( (r<R_crit200[k]*r_limit)&(v<vlimit)&(v>-vlimit) )
		r,v,mags = r[limit],v[limit],mags[limit]
		length.append(len(r))
		### Build Systems ###
		# Method 0
		# Build Ensemble with extra gals to counter shiftgapper
		enc_r0.extend(r[0:gal_num*1.5])
		enc_v0.extend(v[0:gal_num*1.5])
		enc_mag0.extend(mags[0:gal_num*1.5])
		# Build LOS by iterating shiftgapper
		select = []		#Selection Range
		N = 0			#Number of gals
		top = 0			#Highest Index in Select
		go = True		#Loop condition
		while go == True:
			select.append(top)
			top += 1	
	
		line_r0.append(r[0:gal_num])
		line_v0.append(v[0:gal_num])
		line_mag0.append(mags[0:gal_num])




		# Method 1
		rando = random.sample(xrange(len(r)),gal_num)	
				






