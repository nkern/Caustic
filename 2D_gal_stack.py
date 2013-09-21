## This Program uses MDB galaxy data to stack halos and produce a mass estimate

######
from numpy import *
from stack_class_2D import *
from flux_caustics_nideal import *
from matplotlib.pyplot import *
import sys

## DEFINE CONSTANTS ##
h = 0.72 				# Hubble Constant / 100.0
r_limit = 2				# Radius Limit of data in R_crit200 Scale Factor
vlimit = 3500.0				# Velocity Disp limit in km/s
H0 = h*100.0				# Hubble constant
q = 10.0				# Normalization...
c = 300000.0				# Speed of Light km/s
beta = 0.2				# If using constant a beta 

halo_num = 100				# Total number of halos loaded
gal_num = int(sys.argv[2])		# Number of galaxies / halo 
line_num = int(sys.argv[3])		# Number of halo LOS's / ensemble 
bin_num = halo_num/line_num		# Number of bins

run_loc = 'gs_run'+str(sys.argv[4])+''	# Location for data files
arg = [int(sys.argv[1])]		# Ensemble Number (for pbs flux job submission)

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

## Ensemble Arrays
# Diaferio Caustic Mass and Inflection NFW Mass, [Ensemble][Mass]
ENC_DIAMASS = []	
ENC_INFMASS = []
# Caustic and NFW profiles, [Ensemble][ProfileData]
ENC_DIA_CAU = []
ENC_INF_CAU = []
ENC_INF_NFW = []
# R, V and MAG data (stacked), [Ensemble][Data]
ENC_R = []
ENC_V = []
ENC_MAG = []
ENC_GPX3D = []
ENC_GPY3D = []
ENC_GPZ3D = []
ENC_GVX3D = []
ENC_GVY3D = []
ENC_GVZ3D = []
# Ensemble Averaged Properties [Ensemble][Data]
ENC_R200 = []
ENC_M200 = []
ENC_SRAD = []
ENC_ESRAD = []
# Ensemble Calculated Properties [Ensemble][Data]
ENC_VDISP = []

## Individual Line of Sight Arrays
# Line of sight Masses [Ensemble][Line of Sight Data]
LINE_DIAMASS = []
LINE_INFMASS = []
LINE_VDISP = []
LINE_DISSECT = []
# [Ensemble][LOS][Profile]	
LINE_DIACAU = []
##

# Loop over Ensembles
j = 0
for k in arg:

	ENC_R,ENC_V,ENC_MAG,ENC_VDISP,ENC_R200,ENC_M200,ENC_SRAD,ENC_ESRAD,ENC_GPX3D,ENC_GPY3D,ENC_GPZ3D,ENC_GVX3D,ENC_GVY3D,ENC_GVZ3D,LINE_VDISP,LINE_DIAMASS,LINE_INFMASS,LINE_DIACAU,LINE_DISSECT = G.bin_clusters(ENC_R,ENC_V,ENC_MAG,ENC_VDISP,ENC_R200,ENC_M200,ENC_SRAD,ENC_ESRAD,ENC_GPX3D,ENC_GPY3D,ENC_GPZ3D,ENC_GVX3D,ENC_GVY3D,ENC_GVZ3D,LINE_VDISP,Gal_P,Gal_V,Halo_P,Halo_V,M_crit200,R_crit200,SRAD,ESRAD,MAGS,k,r_limit,vlimit,gal_num,line_num,H0,q,c,LINE_DIAMASS,LINE_INFMASS,LINE_DIACAU,LINE_DISSECT,root,beta,scale_data)

	x_range,ENC_DIAMASS,ENC_INFMASS,ENC_DIA_CAU,ENC_INF_CAU,ENC_INF_NFW = G.self_stack_kernel_caustic_masscalc(ENC_R[j],ENC_V[j],ENC_DIAMASS,ENC_INFMASS,ENC_DIA_CAU,ENC_INF_CAU,ENC_INF_NFW,ENC_R200[j],ENC_M200[j],ENC_SRAD[j],ENC_ESRAD[j],ENC_VDISP[j],r_limit,vlimit,H0,q,k,root,beta)

	j += 1	

## Arrays ##
ENC_DIAMASS = array(ENC_DIAMASS)
ENC_INFMASS = array(ENC_INFMASS)
ENC_VDISP = array(ENC_VDISP)
ENC_R,ENC_V,ENC_MAG = array(ENC_R),array(ENC_V),array(ENC_MAG)
ENC_R200,ENC_M200,ENC_SRAD,ENC_ESRAD = array(ENC_R200),array(ENC_M200),array(ENC_SRAD),array(ENC_ESRAD)
ENC_GPX3D,ENC_GPY3D,ENC_GPZ3D,ENC_GVX3D,ENC_GVY3D,ENC_GVZ3D = array(ENC_GPX3D),array(ENC_GPY3D),array(ENC_GPZ3D),array(ENC_GVX3D),array(ENC_GVY3D),array(ENC_GVZ3D)

ENC_DIA_CAU = array(ENC_DIA_CAU)
ENC_INF_CAU = array(ENC_INF_CAU)
ENC_INF_NFW = array(ENC_INF_NFW)

LINE_DIACAU = array(LINE_DIACAU)
LINE_VDISP = array(LINE_VDISP)
LINE_DIAMASS,LINE_INFMASS = array(LINE_DIAMASS),array(LINE_INFMASS)

LINE_DISSECTt = array(LINE_DISSECT)
LINE_DISSECT = zeros(line_num)
for k in range(line_num):
	if k == line_num-1: LINE_DISSECT[k] = -inf	# Just so it will match other line_num length arrays, writing is easier, delete element later
	else: LINE_DISSECT[k] = sum(LINE_DISSECTt[0][:k+1])
#############


raise NameError

if write_data == True:
	## Writing Data ##

	if arg[0] == 0:		# Only write this file once, when it does 0th Ensemble 
		#Writing Constant Data
		f = open(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/program_constants.tab','w')
		f.write(str('#Contains Program Constants'+'\n'))
		f.write(str('#,h,gal_num,line_num,halo_num,bin_num,r_limit,vlimit,beta')+'\n')
		f.write(str(h)+'\t'+str(gal_num)+'\t'+str(line_num)+'\t'+str(halo_num)+'\t'+str(bin_num)+'\t'+str(r_limit)+'\t'+str(vlimit)+'\t'+str(beta)+'\n')
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
	for m in arg:
		#Writing data file 1, containing length 1 arrays
		f = open(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/ensemble_'+str(m)+'_constants.tab','w')
		f.write(str('#Constant data for 2D_gal_self_stack.py run for halo #'+str(m)+'')+'\t'+'\n')
		f.write(str('#')+str('ENC_DIAMASS,ENC_INFMASS,ENC_VDISP,ENC_R200,ENC_M200,ENC_SRAD,ENC_ESRAD,gal_num,line_num')+'\n')
		f.write(str(ENC_DIAMASS[n])+'\t'+str(ENC_INFMASS[n])+'\t'+str(ENC_VDISP[n])+'\t'+str(ENC_R200[n])+'\t'+str(ENC_M200[n])+'\t'+str(ENC_SRAD[n])+'\t'+str(ENC_ESRAD[n])+'\t'+str(gal_num)+'\t'+str(line_num)+'\n')
		f.close()

		#Writing data file 2, containing length line_num arrays
		f = open(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/ensemble_'+str(m)+'_linenum.tab','w')
		f.write(str('#Line of Sigt data for 2D_gal_self_stack.py run for halo #'+str(m)+'')+'\t'+'\n')
		f.write(str('#LINE_DIAMASS,LINE_INFMASS,LINE_VDISP,LINE_DISSECT_NUM')+'\n')
		for j in range(line_num):
			f.write(str(LINE_DIAMASS[n][j])+'\t'+str(LINE_INFMASS[n][j])+'\t'+str(LINE_VDISP[n][j])+'\t'+str(LINE_DISSECT[n][j])+'\n')
		f.close()	

		#Writing data file 3, containing length 201 arrays
		f = open(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/ensemble_'+str(m)+'_profiles.tab','w')
		f.write(str('#Profile data for 2D_gal_self_stack.py run for halo #'+str(m)+'')+'\t'+'\n')
		f.write(str('#')+str('ENC_DIA_CAU')+'\t'+str('ENC_INF_CAU')+'\t'+str('ENC_INF_NFW')+'\t'+str('x_range')+'\n')
		for j in range(201):
			f.write(str(ENC_DIA_CAU[n][j])+'\t'+str(ENC_INF_CAU[n][j])+'\t'+str(ENC_INF_NFW[n][j])+'\t'+str(x_range[j])+'\n')
		f.close()

		#Writing data file 4, containing length line_num*gal_num arrays
		f = open(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/ensemble_'+str(m)+'_RVdata.tab','w')		
		f.write(str('#')+str('Rvalues')+'\t'+str('Vvalues')+'\t'+str('MAGvalues')+'\t'+str('GPX3D')+'\t'+str('GPY3D')+'\t'+str('GPZ3D')+'\t'+str('GVX3D')+'\t'+str('GVY3D')+'\t'+str('GVZ3D')+'\n')
		for j in range(gal_num*line_num):
			f.write(str(ENC_R[n][j])+'\t'+str(ENC_V[n][j])+'\t'+str(ENC_MAG[n][j])+'\t'+str(ENC_GPX3D[n][j])+'\t'+str(ENC_GPY3D[n][j])+'\t'+str(ENC_GPZ3D[n][j])+'\t'+str(ENC_GVX3D[n][j])+'\t'+str(ENC_GVY3D[n][j])+'\t'+str(ENC_GVZ3D[n][j])+'\n')
		f.close()

		#Writing data file 5, containing los profile data
		f = open(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/ensemble_'+str(m)+'_losprofile.tab','w')
		f.write(str('#')+str('Columns [0:line_num] for Dia_Cau,')+'\n')
		for j in range(201):
			for l in range(line_num):
				f.write(str(LINE_DIACAU[n][l][j])+'\t')
			f.write('\n')
		f.close()	
	


	
		n+=1






