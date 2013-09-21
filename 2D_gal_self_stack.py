## This Program uses MDB galaxy data to stack a unique halo multiple times and produce a mass estimate, 2 dimensional

######
print '...importing modules'
from numpy import *
from stack_class_2D import *
from flux_caustics_nideal import *
from matplotlib.pyplot import *
import sys

## DEFINE CONSTANTS ##
h = 1.00 		# Hubble Constant / 100.0
r_limit = 1.25		# Radius Limit of data in R_crit200 Scale Factor
vlimit = 3500.0		# Velocity Disp limit in km/s
H0 = h*100.0		# Hubble constant
q = 10.0		# Normalization...
c = 300000.0		# Speed of Light
halo_num = 100		# Total number of halos loaded
beta = 0.2				# Beta constant
gal_num = int(sys.argv[2])		# Number of galaxies stacked per halo for en. clusters
line_num = int(sys.argv[3])		# Number of lines of sight to stack over
method_num = int(sys.argv[5])		# Which galaxy selection method to use

run_loc = 'ssm'+str(method_num)+'_run'+str(sys.argv[4])+''	# What directory do you want to write files in?

ens_num = [int(sys.argv[1])]	#Ensemble Number, for pbs submission w/ job array

## DEFINE FLAGS ##
use_mems = False	# Membership not known a priori
use_vdisp = True	# If false, do clipping
scale_data = False	# Scale data before stacking

use_flux = True		# Using flux or sophie?
write_data = True	# Write out data?

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


# Diaferio Caustic Mass and Inflection NFW Mass, [Ensemble][Mass]
ENC_CAUMASS = []	
ENC_INFMASS = []
# Caustic and NFW profiles, [Ensemble][ProfileData]
ENC_CAUSURF = []
ENC_INFSURF = []
ENC_INFNFW = []
# R, V and MAG data, [Ensemble][Data]
ENC_R = []
ENC_V = []
ENC_MAG = []
ENC_VDISP = []
ENC_GPX3D = []
ENC_GPY3D = []
ENC_GPZ3D = []
ENC_GVX3D = []
ENC_GVY3D = []
ENC_GVZ3D = []

## Individual Line of Sight Arrays
# Line of sight Masses [Ensemble][Line of Sight Data]
LINE_CAUMASS = []
LINE_INFMASS = []
LINE_VDISP = []
LINE_R = []
LINE_V = []
LINE_MAG = []
# [Ensemble][LOS][Profile]	
LINE_CAUSURF = []
##

# Loop over Ensembles
j = 0
for k in ens_num:

	ENC_R,ENC_V,ENC_MAG,ENC_VDISP,ENC_GPX3D,ENC_GPY3D,ENC_GPZ3D,ENC_GVX3D,ENC_GVY3D,ENC_GVZ3D,LINE_VDISP,LINE_CAUMASS,LINE_INFMASS,LINE_CAUSURF,LINE_R,LINE_V,LINE_MAG = G.self_stack_clusters(ENC_R,ENC_V,ENC_MAG,ENC_VDISP,ENC_GPX3D,ENC_GPY3D,ENC_GPZ3D,ENC_GVX3D,ENC_GVY3D,ENC_GVZ3D,LINE_VDISP,Gal_P,Gal_V,Halo_P,Halo_V,M_crit200,R_crit200,SRAD,ESRAD,MAGS,k,r_limit,vlimit,gal_num,line_num,method_num,H0,q,c,LINE_CAUMASS,LINE_INFMASS,LINE_CAUSURF,LINE_R,LINE_V,LINE_MAG,root,beta)

	x_range,ENC_CAUMASS,ENC_INFMASS,ENC_CAUSURF,ENC_INFSURF,ENC_INFNFW = G.self_stack_kernel_caustic_masscalc(ENC_R[j],ENC_V[j],ENC_CAUMASS,ENC_INFMASS,ENC_CAUSURF,ENC_INFSURF,ENC_INFNFW,R_crit200[k],M_crit200[k],SRAD[k],ESRAD[k],ENC_VDISP[j],r_limit,vlimit,H0,q,k,root,beta)

	j += 1	

## Arrays ##
ENC_CAUMASS = array(ENC_CAUMASS)
ENC_INFMASS = array(ENC_INFMASS)
ENC_VDISP = array(ENC_VDISP)
ENC_R,ENC_V,ENC_MAG = array(ENC_R),array(ENC_V),array(ENC_MAG)
ENC_GPX3D,ENC_GPY3D,ENC_GPZ3D,ENC_GVX3D,ENC_GVY3D,ENC_GVZ3D = array(ENC_GPX3D),array(ENC_GPY3D),array(ENC_GPZ3D),array(ENC_GVX3D),array(ENC_GVY3D),array(ENC_GVZ3D)

ENC_CAUSURF = array(ENC_CAUSURF)
ENC_INFSURF = array(ENC_INFSURF)
ENC_INFNFW = array(ENC_INFNFW)

LINE_CAUSURF = array(LINE_CAUSURF)
LINE_VDISP = array(LINE_VDISP)
LINE_CAUMASS,LINE_INFMASS = array(LINE_CAUMASS),array(LINE_INFMASS)

LINE_R,LINE_V,LINE_MAG = array(LINE_R),array(LINE_V),array(LINE_MAG)
#############

if write_data == True:
	## Writing Data ##

	if ens_num[0] == 0:
		#Writing Constant Data
		f = open(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/program_constants.tab','w')
		f.write(str('#Contains Program Constants'+'\n'))
		f.write(str('#')+str('h')+'\t'+str('gal_num')+'\t'+str('line_num')+'\t'+str('halo_num')+'\t'+str('method_num')+'\t'+str('r_limit')+'\t'+str('vlimit')+'\t'+str('beta')+'\n')
		f.write(str(h)+'\t'+str(gal_num)+'\t'+str(line_num)+'\t'+str(halo_num)+'\t'+str(method_num)+'\t'+str(r_limit)+'\t'+str(vlimit)+'\t'+str(beta)+'\n')
		f.close()

		#Writing More Constant Data, containing halo_num length arrays
		f = open(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/simdata.tab','w')
		f.write(str('#Contains Table Value Halo Data Arrays')+'\n')
		f.write(str('#HaloID,M_crit200,R_crit200,SRAD,ESRAD,HVD')+'\n')
		for j in range(halo_num):
			f.write(str(HaloID[j])+'\t'+str(M_crit200[j])+'\t'+str(R_crit200[j])+'\t'+str(SRAD[j])+'\t'+str(ESRAD[j])+'\t'+str(HVD[j])+'\n')
		f.close()

	#Writing Halo Data
	n = 0
	for m in ens_num:
		#Writing data file 1, containing length 1 arrays
		f = open(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/halo_'+str(m)+'_constants.tab','w')
		f.write(str('#Constant data for 2D_gal_self_stack.py run for halo #'+str(m)+'')+'\t'+'\n')
		f.write(str('#')+str('ENC_CAUMASS,ENC_INFMASS,ENC_VDISP,gal_num,line_num,method_num')+'\n')
		f.write(str(ENC_CAUMASS[n])+'\t'+str(ENC_INFMASS[n])+'\t'+str(ENC_VDISP[n])+'\t'+str(gal_num)+'\t'+str(line_num)+'\t'+str(method_num)+'\n')
		f.close()

		#Writing data file 2, containing length line_num arrays
		f = open(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/halo_'+str(m)+'_linenum.tab','w')
		f.write(str('#Line of Sigt data for 2D_gal_self_stack.py run for halo #'+str(m)+'')+'\t'+'\n')
		f.write(str('#')+str('LINE_CAUMASS,LINE_INFMASS,LINE_VDISP')+'\n')
		for j in range(line_num):
			f.write(str(LINE_CAUMASS[n][j])+'\t'+str(LINE_INFMASS[n][j])+'\t'+str(LINE_VDISP[n][j])+'\n')
		f.close()	

		#Writing data file 3, containing length 201 arrays
		f = open(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/halo_'+str(m)+'_profiles.tab','w')
		f.write(str('#Profile data for 2D_gal_self_stack.py run for halo #'+str(m)+'')+'\t'+'\n')
		f.write(str('#')+str('ENC_CAUSURF')+'\t'+str('ENC_INFSURF')+'\t'+str('ENC_INFNFW')+'\t'+str('x_range')+'\n')
		for j in range(201):
			f.write(str(ENC_CAUSURF[n][j])+'\t'+str(ENC_INFSURF[n][j])+'\t'+str(ENC_INFNFW[n][j])+'\t'+str(x_range[j])+'\n')
		f.close()

		#Writing data file 4, containing length line_num*gal_num arrays
		f = open(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/halo_'+str(m)+'_RVdata.tab','w')		
		f.write(str('#')+str('Rvalues')+'\t'+str('Vvalues')+'\t'+str('MAGvalues')+'\t'+str('GPX3D')+'\t'+str('GPY3D')+'\t'+str('GPZ3D')+'\t'+str('GVX3D')+'\t'+str('GVY3D')+'\t'+str('GVZ3D')+'\n')
		for j in range( len(ENC_R[n]) ):
			f.write(str(ENC_R[n][j])+'\t'+str(ENC_V[n][j])+'\t'+str(ENC_MAG[n][j])+'\t'+str(ENC_GPX3D[n][j])+'\t'+str(ENC_GPY3D[n][j])+'\t'+str(ENC_GPZ3D[n][j])+'\t'+str(ENC_GVX3D[n][j])+'\t'+str(ENC_GVY3D[n][j])+'\t'+str(ENC_GVZ3D[n][j])+'\n')
		f.close()

		#Writing data file 5, containing los profile data
		f = open(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/halo_'+str(m)+'_losprofile.tab','w')
		f.write(str('#')+str('Columns [0:line_num] for Caustic surface,')+'\n')
		for j in range(201):
			for l in range(line_num):
				f.write(str(LINE_CAUSURF[n][l][j])+'\t')
			f.write('\n')
		f.close()	
	
		#Writing data file(s) 6, containing LOS RV data (arrays vary in length)
		for l in range(line_num):
			f = open(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+run_loc+'/LOS_RV/halo_'+str(m)+'_los_'+str(l)+'_rv.tab','w')
			f.write(str('#Contains R,V,MAG data for a particular LOS')+'\n')
			f.write(str('#line_r,line_v,line_mag')+'\n')
			for j in range( LINE_R[n][l].shape[0]):
				f.write( str(LINE_R[n][l][j])+'\t'+str(LINE_V[n][l][j])+'\t'+str(LINE_MAG[n][l][j])+'\n')
			f.close()	

		n+=1


