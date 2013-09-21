## This program takes the stacking code and loops them to do different things

#Program 1: Calculates Bias as a function of stacked halos using galaxy data
#Program 2: Looks at Phase Space plots for individiual halos and ensemble clusters for GALAXIES
#Program 3: Looks at Phase Space plots for individual halos and ensemble clusters for PARTICLES
#Program 4: Loops Program 3 to look at how unique runs (particles) produce unique results for caustic surface and mass estimates, and attempts to quantify it via standard dev. 

## IMPORT MODULES ##

import os
from numpy import *
from stack_class_3D import stack_3D
from matplotlib.pyplot import *

## DEFINE FLAGS ##

run_code1 = False		# Which program are you running? Only 1 should == True.
run_code2 = True
run_code3 = False
run_code4 = False

## DEFINE CONSTANTS ##


## INITIALIZATION ##

Calc = stack_3D()

## USER FUNCTIONS ##

def clearall():
	"""clear all globals, actually not necessary when using a function"""
	for uniquevar in [var for var in globals().copy() if var[0] != '_' and var != 'clearall' and var != 'arange' and var != 'BIAS1' and var != 'BIAS2' and var != 'BIAS3' and var != 'bin_range2']:
		del globals()[uniquevar]
	return

##################



## PROGRAM 1 ##
if run_code1 == True:
	print ""
	print "####################################"
	print "##### RUNNING 3D_calc_stack.py #####"
	print "####################################"
	print ""

	## DEFINE CONSTANTS AND FLAGS ##
	gal_num = 1000		# number of galaxies per halo taken when stacked
	code_num = 1 		# program number, used in assigning return values in stack_class_3D.py
	run_num = [0,100]
	##
	BIAS1 = []		# different ways to calculate bias
	BIAS2 = []
	BIAS3 = []		
	bin_range2 = []		# holds the bin_range value
	j = 0
	for k in [2,4,5,10,20,25,50]:	# All multiples of 100 so that every halo is used in the calculations

		print "############ LOOP "+str(j)+", bin_range = "+str(k)+" ############"
		bias1,bias2,bias3 = Calc.gal_stack_3D(k,gal_num,run_num,code_num)
		BIAS1.append(bias1)
		BIAS2.append(bias2)
		BIAS3.append(bias3)
		bin_range2.append(k)	
		j += 1

	BIAS1,BIAS2,BIAS3,bin_range2 = array(BIAS1),array(BIAS2),array(BIAS3),array(bin_range2)

	print ""
	print "## DONE WITH 3D_calc_stack.py ##"

###########################



## PROGRAM 2 ##
if run_code2 == True:
	print ""
	print "####################################"
	print "##### RUNNING 3D_calc_stack.py #####"
	print "#####      Using Galaxies      #####"
	print "####################################"
	print ""

	## DEFINE FLAGS AND CONSTANTS ##

	use_gals = True		# using Galaxies!!... outdated..	
	code_num = 2		# Program number, used in assigning return values in stack_class_3D.py

	##########
	print ""
	print "## WORKING ON INDIVIDUAL HALOS ##"
	print ""	
	## Individual Halos ##
	gal_num = 100		# number of galaxies per halo to take
	run_num_ha = [0,100]	# run number for halos

	x_range,INF_NFWMASS,DIA_NFWMASS,INF_CAUMASS,DIA_CAUMASS,INF_MPROF,INF_NFW,INF_CAU,DIA_MPROF,DIA_NFW,DIA_CAU,R,V,MAGS,M_crit200,R_crit200 = Calc.halo_gal_3D(gal_num,run_num_ha,code_num)		
	##########
	print ""
	print "## WORKING ON ENSEMBLE CLUSTERS ##"
	print ""
	## Ensemble Clusters ##
	bin_range = 10		# Number of halos per ensemble cluster
	gal_num = 10		# Number of galaxies per halo to take
	run_num_en = [0,10]	# run number for ensembles

	ENC_x_range,ENC_INF_NFWMASS,ENC_DIA_NFWMASS,ENC_INF_CAUMASS,ENC_DIA_CAUMASS,ENC_INF_MPROF,ENC_INF_NFW,ENC_INF_CAU,ENC_DIA_MPROF,ENC_DIA_NFW,ENC_DIA_CAU,ENC_M200,ENC_R200,ENC_R,ENC_V,ENC_MAG = Calc.gal_stack_3D(bin_range,gal_num,run_num_en,code_num)



	
	
## PROGRAM 3 ##
if run_code3 == True:		

	print ""
	print "####################################"
	print "##### RUNNING 3D_calc_stack.py #####"
	print "#####	Using Particles	      #####"
	print "####################################"
	print ""

	## DEFINE FLAGS AND CONSTANTS ##

	use_gals = False	# using Particles!!	
	code_num = 3		# Program number, used in assigning return values in stack_class_3D.py

	##########
	print ""
	print "## WORKING ON INDIVIDUAL HALOS ##"
	print ""	
	## Individual Halos ##
	gal_num = 10000				# number of galaxies per halo to take
	run_num = [0,100]			# run number for halos

	x_range,INF_NFWMASS,DIA_NFWMASS,INF_CAUMASS,DIA_CAUMASS,INF_MPROF,INF_NFW,INF_CAU,DIA_MPROF,DIA_NFW,DIA_CAU,R,V,HaloID,R_crit200,M_crit200 = Calc.halo_part_3D(gal_num,run_num,code_num)	


	##########
	print ""
	print "## WORKING ON ENSEMBLE CLUSTERS ##"
	print ""
	## Ensemble Clusters ##
	bin_range = 10				# Number of halos per ensemble cluster
	gal_num = 1000				# Number of particles per halo to take

	ENC_x_range,ENC_INF_NFWMASS,ENC_DIA_NFWMASS,ENC_INF_CAUMASS,ENC_DIA_CAUMASS,ENC_INF_MPROF,ENC_INF_NFW,ENC_INF_CAU,ENC_DIA_MPROF,ENC_DIA_NFW,ENC_DIA_CAU,ENC_R,ENC_V,ENC_M200,ENC_R200 = Calc.part_stack_3D(bin_range,gal_num,run_num,code_num)

	## creating averaged arrays ##
	ENC_INF_CAUMEAN = []
	ENC_INF_NFWMEAN = []
	ENC_DIA_CAUMEAN = []
	ENC_DIA_NFWMEAN = []

	inf_nfw = copy(INF_NFW)
	inf_cau = copy(INF_CAU)
	dia_nfw = copy(DIA_NFW)
	dia_cau = copy(DIA_CAU)
	for j in range( (run_num[1]-run_num[0])/bin_range ):
		inf_caumean = []
		inf_nfwmean = []
		dia_caumean = []
		dia_nfwmean = []

		for k in range(201):
			inf_caumean.append( mean( inf_cau[j*bin_range:(j+1)*bin_range].T[k] ) )
			inf_nfwmean.append( mean( inf_nfw[j*bin_range:(j+1)*bin_range].T[k] ) )
			dia_caumean.append( mean( dia_cau[j*bin_range:(j+1)*bin_range].T[k] ) )
			dia_nfwmean.append( mean( dia_nfw[j*bin_range:(j+1)*bin_range].T[k] ) )  
	
		ENC_INF_CAUMEAN.append( array(inf_caumean)  )
		ENC_INF_NFWMEAN.append( array(inf_nfwmean) )
		ENC_DIA_CAUMEAN.append( array(dia_caumean) )
		ENC_DIA_NFWMEAN.append( array(dia_nfwmean) )

	ENC_INF_CAUMEAN,ENC_INF_NFWMEAN = array(ENC_INF_CAUMEAN),array(ENC_INF_NFWMEAN)
	ENC_DIA_CAUMEAN,ENC_DIA_NFWMEAN = array(ENC_DIA_CAUMEAN),array(ENC_DIA_NFWMEAN)

#############################



## PROGRAM 4 ##

if run_code4 == True:		

	tot_enc_inf_cau = []
	tot_enc_inf_nfw = []
	tot_enc_dia_cau = []
	tot_enc_dia_nfw = []

	tot_inf_cau = []
	tot_inf_nfw = []
	tot_dia_cau = []
	tot_dia_nfw = []

	tot_enc_inf_caumean = []
	tot_enc_inf_nfwmean = []
	tot_enc_dia_caumean = []
	tot_enc_dia_nfwmean = []

	tot_inf_nfwmass = []	
	tot_dia_caumass = []
	tot_enc_inf_nfwmass = []
	tot_enc_dia_caumass = []

	for j in range(5):
		

		print ""
		print "####################################"
		print "##### RUNNING 3D_calc_stack.py #####"
		print "#####	Using Particles	      #####"
		print "####################################"
		print ""

		## DEFINE FLAGS AND CONSTANTS ##

		use_gals = False	# using Particles!!	
		code_num = 4		# Program number, used in assigning return values in stack_class_3D.py

		##########
		print ""
		print "## WORKING ON INDIVIDUAL HALOS ##"
		print ""	
		## Individual Halos ##
		gal_num = 10000				# number of galaxies per halo to take
		run_num = [0,10]			# run number for halos

		x_range,INF_NFWMASS,DIA_NFWMASS,INF_CAUMASS,DIA_CAUMASS,INF_MPROF,INF_NFW,INF_CAU,DIA_MPROF,DIA_NFW,DIA_CAU,R,V,HaloID,R_crit200,M_crit200 = Calc.halo_part_3D(gal_num,run_num,code_num)	


		##########
		print ""
		print "## WORKING ON ENSEMBLE CLUSTERS ##"
		print ""
		## Ensemble Clusters ##
		bin_range = 10				# Number of halos per ensemble cluster
		gal_num = 1000				# Number of particles per halo to take

		ENC_x_range,ENC_INF_NFWMASS,ENC_DIA_NFWMASS,ENC_INF_CAUMASS,ENC_DIA_CAUMASS,ENC_INF_MPROF,ENC_INF_NFW,ENC_INF_CAU,ENC_DIA_MPROF,ENC_DIA_NFW,ENC_DIA_CAU,ENC_R,ENC_V,ENC_M200,ENC_R200 = Calc.part_stack_3D(bin_range,gal_num,run_num,code_num)

		## creating averaged arrays ##
		ENC_INF_CAUMEAN = []
		ENC_INF_NFWMEAN = []
		ENC_DIA_CAUMEAN = []
		ENC_DIA_NFWMEAN = []

		inf_nfw = copy(INF_NFW)
		inf_cau = copy(INF_CAU)
		dia_nfw = copy(DIA_NFW)
		dia_cau = copy(DIA_CAU)
		for j in range( (run_num[1]-run_num[0])/bin_range ):
			inf_caumean = []
			inf_nfwmean = []
			dia_caumean = []
			dia_nfwmean = []

			for k in range(201):
				inf_caumean.append( mean( inf_cau[j*bin_range:(j+1)*bin_range].T[k] ) )
				inf_nfwmean.append( mean( inf_nfw[j*bin_range:(j+1)*bin_range].T[k] ) )
				dia_caumean.append( mean( dia_cau[j*bin_range:(j+1)*bin_range].T[k] ) )
				dia_nfwmean.append( mean( dia_nfw[j*bin_range:(j+1)*bin_range].T[k] ) )  
		
			ENC_INF_CAUMEAN.append( array(inf_caumean)  )
			ENC_INF_NFWMEAN.append( array(inf_nfwmean) )
			ENC_DIA_CAUMEAN.append( array(dia_caumean) )
			ENC_DIA_NFWMEAN.append( array(dia_nfwmean) )

		ENC_INF_CAUMEAN,ENC_INF_NFWMEAN = array(ENC_INF_CAUMEAN),array(ENC_INF_NFWMEAN)
		ENC_DIA_CAUMEAN,ENC_DIA_NFWMEAN = array(ENC_DIA_CAUMEAN),array(ENC_DIA_NFWMEAN)

		tot_enc_inf_cau.append(ENC_INF_CAU)
		tot_enc_inf_nfw.append(ENC_INF_NFW)
		tot_enc_dia_cau.append(ENC_DIA_CAU)
		tot_enc_dia_nfw.append(ENC_DIA_NFW)
	
		tot_inf_cau.append(INF_CAU)
		tot_inf_nfw.append(INF_NFW)
		tot_dia_cau.append(DIA_CAU)
		tot_dia_nfw.append(DIA_NFW)

		tot_enc_inf_caumean.append(ENC_INF_CAUMEAN)
		tot_enc_inf_nfwmean.append(ENC_INF_NFWMEAN)
		tot_enc_dia_caumean.append(ENC_DIA_CAUMEAN)
		tot_enc_dia_nfwmean.append(ENC_DIA_NFWMEAN)

		tot_inf_nfwmass.append(INF_NFWMASS)
		tot_dia_caumass.append(DIA_CAUMASS)
		tot_enc_inf_nfwmass.append(ENC_INF_NFWMASS)
		tot_enc_dia_caumass.append(ENC_DIA_CAUMASS)
	



