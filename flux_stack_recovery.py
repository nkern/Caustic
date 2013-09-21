'''This program takes files outputted by 2D_gal_self_stack.py, and loads in the data to interpreter to manipulate'''

from numpy import *
from matplotlib.pyplot import *
import matplotlib.mlab as mlab
import numpy.ma as ma
import astStats

## FLAGS ##
use_flux = True 	# Using flux or sophie?
get_los = False		# Upload los rv data (takes a while)


run_loc = 'self_stack_run2'	#Which data directory do you want to open, *run1 or *run2?

if use_flux == True:
	root = str('/nfs/christoq_ls')
else:
	root = str('/n/Christoq1')


## FUNCTIONS ##

def ss_recover():
	#Preliminary data file upload
	global h,gal_num,line_num,halo_num,r_limit,vlimit,beta
	h,gal_num,line_num,halo_num,r_limit,vlimit,beta = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/program_constants.tab',unpack=True)
	halo_num = int(halo_num)
	line_num,gal_num = int(line_num),int(gal_num)

	#Second preliminary data file upload
	global HaloID,M_crit200,R_crit200,SRAD,ESARD,HVD
	HaloID,M_crit200,R_crit200,SRAD,ESRAD,HVD = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/simdata.tab',unpack=True)	
	HaloID = str(HaloID)
	HaloID,M_crit200,R_crit200,SRAD,ESRAD,HVD = HaloID[:halo_num],M_crit200[:halo_num],R_crit200[:halo_num],SRAD[:halo_num],ESRAD[:halo_num],HVD[:halo_num]


	#First Data file upload
	global ENC_CAUMASS,ENC_INFMASS,ENC_VDISP
	j = 0
	for m in range(halo_num):
		if j == 0:	#Initialization of arrays
			ENC_CAUMASS,ENC_INFMASS,ENC_VDISP = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/halo_'+str(m)+'_constants.tab',usecols=(0,1,2),unpack=True)
		else:
			ENC_CAUMASSt,ENC_INFMASSt,ENC_VDISPt = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/halo_'+str(m)+'_constants.tab',usecols=(0,1,2),unpack=True)
			ENC_CAUMASS = hstack([ENC_CAUMASS,ENC_CAUMASSt])
			ENC_INFMASS = hstack([ENC_INFMASS,ENC_INFMASSt])
			ENC_VDISP = hstack([ENC_VDISP,ENC_VDISPt])
		j +=1

	#Second data file upload
	global LINE_CAUMASS,LINE_INFMASS,LINE_VDISP
	j = 0
	for m in range(halo_num):
		if j == 0:	#Initialization of arrays
			LINE_CAUMASS,LINE_INFMASS,LINE_VDISP = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/halo_'+str(m)+'_linenum.tab',unpack=True)
		else:
			line_caumass,line_infmass,line_vdisp = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/halo_'+str(m)+'_linenum.tab',unpack=True)
			LINE_CAUMASS = vstack([LINE_CAUMASS,line_caumass])
			LINE_INFMASS = vstack([LINE_INFMASS,line_infmass])
			LINE_VDISP = vstack([LINE_VDISP,line_vdisp])
		j+=1

	#Third data file upload
	global ENC_CAUSURF,ENC_INFSURF,ENC_INFNFW,x_range
	j = 0
	for m in range(halo_num):
		if j == 0:	#Initialization of arrays
			ENC_CAUSURF, ENC_INFSURF, ENC_INFNFW, x_range = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/halo_'+str(m)+'_profiles.tab',unpack=True)
		else:
			enc_causurf,enc_infsurf,enc_infnfw, x_range = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/halo_'+str(m)+'_profiles.tab',unpack=True)
			ENC_CAUSURF = vstack([ENC_CAUSURF,enc_causurf])
			ENC_INFSURF = vstack([ENC_INFSURF,enc_infsurf])
			ENC_INFNFW = vstack([ENC_INFNFW,enc_infnfw])
		j+=1

	#Fourth data file upload
	global ENC_R,ENC_V,ENC_MAG,ENC_GPX3D,ENC_GPY3D,ENC_GPZ3D,ENC_GVX3D,ENC_GVY3D,ENC_GVZ3D
	ENC_R,ENC_V,ENC_MAG,ENC_GPX3D,ENC_GPY3D,ENC_GPZ3D,ENC_GVX3D,ENC_GVY3D,ENC_GVZ3D = [],[],[],[],[],[],[],[],[]
	j = 0
	for m in range(halo_num):
		enc_r,enc_v,enc_mag,enc_gpx3d,enc_gpy3d,enc_gpz3d,enc_gvx3d,enc_gvy3d,enc_gvz3d = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/halo_'+str(m)+'_RVdata.tab',unpack=True)
		ENC_R.append(enc_r)
		ENC_V.append(enc_v) 
		ENC_MAG.append(enc_mag) 
		ENC_GPX3D.append(enc_gpx3d) 
		ENC_GPY3D.append(enc_gpy3d) 
		ENC_GPZ3D.append(enc_gpz3d)
		ENC_GVX3D.append(enc_gvx3d) 
		ENC_GVY3D.append(enc_gvy3d) 
		ENC_GVZ3D.append(enc_gvz3d) 
		j += 1
	ENC_R,ENC_V,ENC_MAG,ENC_GPX3D,ENC_GPY3D,ENC_GPZ3D,ENC_GVX3D,ENC_GVY3D,ENC_GVZ3D = array(ENC_R),array(ENC_V),array(ENC_MAG),array(ENC_GPX3D),array(ENC_GPY3D),array(ENC_GPZ3D),array(ENC_GVX3D),array(ENC_GVY3D),array(ENC_GVZ3D)  	
	#Fifth data file to upload
	global LINE_CAUSURF
	j = 0
	for m in range(halo_num):
		if j == 0:
			line_prof = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/halo_'+str(m)+'_losprofile.tab',unpack=True)
			LINE_CAUSURF = array([line_prof[0:line_num]])
		else:
			line_prof = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/halo_'+str(m)+'_losprofile.tab',unpack=True)
			line_causurf = array([line_prof[0:line_num]])
			LINE_CAUSURF = vstack([LINE_CAUSURF,line_causurf])
		j += 1

	#Sixth data set upload (los rv data)
	if get_los == True:
		global LINE_R, LINE_V, LINE_MAG
		LINE_R,LINE_V,LINE_MAG = [],[],[]
		j = 0
		for m in range(halo_num):	
			line_r,line_v,line_mag = [],[],[]
			for l in range(line_num):	
					r,v,mag = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/LOS_RV/halo_'+str(m)+'_los_'+str(l)+'_rv.tab',unpack=True)
					line_r.append(r)
					line_v.append(v)
					line_mag.append(mag)
			LINE_R.append(line_r)
			LINE_V.append(line_v)
			LINE_MAG.append(line_mag)
		LINE_R,LINE_V,LINE_MAG = array(LINE_R),array(LINE_V),array(LINE_MAG)

	# Other data arrays to use:
	global avg_mfrac,avg_hvdfrac,stack_mfrac,stack_hvdfrac,maLINE_CAUMASS,maLINE_VDISP
	global stack_mbias,stack_mscat,stack_vbias,stack_vscat,avg_mbias,avg_mscat,avg_vbias,avg_vscat

	maLINE_CAUMASS = ma.masked_array(LINE_CAUMASS,mask=LINE_CAUMASS==0)	#Mask 0 Values
	maLINE_VDISP = ma.masked_array(LINE_VDISP,mask=LINE_VDISP==0)		#Mask 0 Values

	### Mass Fractions ###
	# Note: I was using map() as an iterator, but for N = 5, sometimes there are less than 3 non-masked values per los
	# Note: and biweight###() does not take less than 4 unique values. I don't yet know how to incorporate a "try:" 
	# Note: statement into an iterator function like map(), so I resort to a "for" loop
	## Ensemble fractions
	stack_mfrac = ma.log(ENC_CAUMASS/M_crit200)
	stack_hvdfrac = ma.log(ENC_VDISP/HVD)
	## Averaged fractions
	a_size = halo_num		# This becomes line_num if doing vertical average first!!
	avg_mfrac,avg_hvdfrac = zeros(a_size),zeros(a_size)
	for a in range(a_size):	
		try:
			avg_mfrac[a] = astStats.biweightLocation( ma.copy(ma.log(maLINE_CAUMASS[a]/M_crit200[a])), 6.0 )
			avg_hvdfrac[a] = astStats.biweightLocation( ma.copy(ma.log(maLINE_VDISP[a]/HVD[a])), 6.0 )
		except:
			avg_mfrac[a] = ma.mean( ma.log(maLINE_CAUMASS[a]/M_crit200[a]) )
			avg_hvdfrac[a] = ma.mean( ma.log(maLINE_VDISP[a]/M_crit200[a]) )	
	# Bias and Scatter for Ensemble and LOS Average Systems
	stack_mbias,stack_mscat = astStats.biweightLocation(ma.copy(stack_mfrac),6.0),astStats.biweightScale(ma.copy(stack_mfrac),9.0)
	avg_mbias,avg_mscat = astStats.biweightLocation(ma.copy(avg_mfrac),6.0),astStats.biweightScale(ma.copy(avg_mfrac),9.0)
	stack_vbias,stack_vscat = astStats.biweightLocation(ma.copy(stack_hvdfrac),6.0),astStats.biweightScale(ma.copy(stack_hvdfrac),9.0)
	avg_vbias,avg_vscat = astStats.biweightLocation(ma.copy(avg_hvdfrac),6.0),astStats.biweightScale(ma.copy(avg_hvdfrac),9.0)	 

		#########################################
	
		#########################################


def gs_recover():
	#Preliminary data file upload
	global h,gal_num,line_num,halo_num,bin_num,r_limit,vlimit,beta
	h,gal_num,line_num,halo_num,bin_num,r_limit,vlimit,beta = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/program_constants.tab',unpack=True)
	halo_num = int(halo_num)
	line_num,gal_num,bin_num = int(line_num),int(gal_num),int(bin_num)

	#Second preliminary data file upload
	global HaloID,M_crit200,R_crit200,SRAD,ESARD,HVD
	HaloID,M_crit200,R_crit200,SRAD,ESRAD,HVD = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/simdata.tab',unpack=True)	
	HaloID = str(HaloID)

	#First Data file upload
	global ENC_DIAMASS,ENC_INFMASS,ENC_VDISP,ENC_R200,ENC_M200,ENC_SRAD,ENC_ESRAD
	j = 0
	for m in range(bin_num):
		if j == 0:	#Initialization of arrays
			ENC_DIAMASS,ENC_INFMASS,ENC_VDISP,ENC_R200,ENC_M200,ENC_SRAD,ENC_ESRAD = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/ensemble_'+str(m)+'_constants.tab',usecols=(0,1,2,3,4,5,6),unpack=True)
		else:
			ENC_DIAMASSt,ENC_INFMASSt,ENC_VDISPt,ENC_R200t,ENC_M200t,ENC_SRADt,ENC_ESRADt = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/ensemble_'+str(m)+'_constants.tab',usecols=(0,1,2,3,4,5,6),unpack=True)
			ENC_DIAMASS = hstack([ENC_DIAMASS,ENC_DIAMASSt])
			ENC_INFMASS = hstack([ENC_INFMASS,ENC_INFMASSt])
			ENC_VDISP = hstack([ENC_VDISP,ENC_VDISPt])
			ENC_R200 = hstack([ENC_R200,ENC_R200t])
			ENC_M200 = hstack([ENC_M200,ENC_M200t])
			ENC_SRAD = hstack([ENC_SRAD,ENC_SRADt])
			ENC_ESRAD = hstack([ENC_ESRAD,ENC_ESRADt])

		j +=1

	#Second data file upload
	global LINE_DIAMASS,LINE_INFMASS,LINE_VDISP,LINE_DISSECT
	j = 0
	for m in range(bin_num):
		if j == 0:	#Initialization of arrays
			LINE_DIAMASS,LINE_INFMASS,LINE_VDISP,LINE_DISSECT = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/ensemble_'+str(m)+'_linenum.tab',unpack=True)
			LINE_DISSECT = delete(LINE_DISSECT,-1)
		else:
			line_diamass,line_infmass,line_vdisp,line_dissect = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/ensemble_'+str(m)+'_linenum.tab',unpack=True)
			line_dissect = delete(line_dissect,-1)
			LINE_DIAMASS = vstack([LINE_DIAMASS,line_diamass])
			LINE_INFMASS = vstack([LINE_INFMASS,line_infmass])
			LINE_VDISP = vstack([LINE_VDISP,line_vdisp])
			LINE_DISSECT = vstack([LINE_DISSECT,line_dissect])	
		j+=1

	#Third data file upload
	global ENC_DIACAU,ENC_INF_CAU,ENC_INF_NFW,x_range
	j = 0
	for m in range(bin_num):
		if j == 0:	#Initialization of arrays
			ENC_DIACAU, ENC_INF_CAU, ENC_INF_NFW, x_range = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/ensemble_'+str(m)+'_profiles.tab',unpack=True)
		else:
			enc_dia_cau,enc_inf_cau,enc_inf_nfw, x_range = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/ensemble_'+str(m)+'_profiles.tab',unpack=True)
			ENC_DIACAU = vstack([ENC_DIACAU,enc_dia_cau])
			ENC_INF_CAU = vstack([ENC_INF_CAU,enc_inf_cau])
			ENC_INF_NFW = vstack([ENC_INF_NFW,enc_inf_nfw])
		j+=1

	#Fourth data file upload
	global ENC_R,ENC_V,ENC_MAG,ENC_GPX3D,ENC_GPY3D,ENC_GPZ3D,ENC_GVX3D,ENC_GVY3D,ENC_GVZ3D,LINE_R,LINE_V,LINE_MAG
	j = 0
	for m in range(bin_num):
		if j == 0: 	#Initialization of arrays
			ENC_R,ENC_V,ENC_MAG,ENC_GPX3D,ENC_GPY3D,ENC_GPZ3D,ENC_GVX3D,ENC_GVY3D,ENC_GVZ3D = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/ensemble_'+str(m)+'_RVdata.tab',unpack=True)
			LINE_R = array([hsplit(ENC_R,LINE_DISSECT[m])])
			LINE_V = array([hsplit(ENC_V,LINE_DISSECT[m])])
			LINE_MAG = array([hsplit(ENC_MAG,LINE_DISSECT[m])])
		else:
			enc_r,enc_v,enc_mag,enc_gpx3d,enc_gpy3d,enc_gpz3d,enc_gvx3d,enc_gvy3d,enc_gvz3d = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/ensemble_'+str(m)+'_RVdata.tab',unpack=True)
			line_r,line_v,line_mag = array([hsplit(enc_r,LINE_DISSECT[m])]),array([hsplit(enc_v,LINE_DISSECT[m])]),array([hsplit(enc_mag,LINE_DISSECT[m])])
			ENC_R = vstack([ENC_R,enc_r])	
			ENC_V = vstack([ENC_V,enc_v])
			ENC_MAG = vstack([ENC_MAG,enc_mag])
			ENC_GPX3D = vstack([ENC_GPX3D,enc_gpx3d])
			ENC_GPY3D = vstack([ENC_GPY3D,enc_gpy3d])
			ENC_GPZ3D = vstack([ENC_GPZ3D,enc_gpz3d])
			ENC_GVX3D = vstack([ENC_GVX3D,enc_gvx3d])
			ENC_GVY3D = vstack([ENC_GVY3D,enc_gvy3d])
			ENC_GVZ3D = vstack([ENC_GVZ3D,enc_gvz3d])
			LINE_R = vstack([LINE_R,line_r])
			LINE_V = vstack([LINE_V,line_v])
			LINE_MAG = vstack([LINE_MAG,line_mag])
		j += 1

	#Fifth data file to upload
	global LINE_DIACAU
	j = 0
	for m in range(bin_num):
		if j == 0:
			line_prof = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/ensemble_'+str(m)+'_losprofile.tab',unpack=True)
			LINE_DIACAU = array([line_prof[0:line_num]])
		else:
			line_prof = loadtxt(''+root+'/nkern/Documents/MDB_milliMil_halodata/Caustic/stack_data/'+str(run_loc)+'/ensemble_'+str(m)+'_losprofile.tab',unpack=True)
			line_diacau = array([line_prof[0:line_num]])
			LINE_DIACAU = vstack([LINE_DIACAU,line_diacau])
		j += 1

	# Other data arrays to use:
	global mean_line_diamass,mean_line_vdisp
	mean_line_diamass = zeros(bin_num)
	mean_line_vdisp = zeros(bin_num)
	for j in range(bin_num):
		mean_line_diamass[j] = mean(LINE_DIAMASS[j])
		mean_line_vdisp[j] = mean(LINE_VDISP[j])


