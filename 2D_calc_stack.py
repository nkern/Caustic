'''This program takes Stacking Code for Projected Space and loops them to calculate different things'''

# Program 1: calculates inherent Caustic Method bias as a function of 'sphericity'
# Program 2: looks at phase space plots and 1-1 lines for halos and ensembles, GALAXIES
# Program 3: looks at phase space plots and 1-1 lines for halos and ensembles, PARTICLES
# Program 4: single halo stacking, or self-stack.

## IMPORT MODULES ##

import os
from numpy import *
from stack_class_2D import *
from matplotlib.pyplot import *

## INITIALIZATION ##
Calc = stack_2D()

## DEFINE FLAGS ##

run_code1 = False
run_code2 = True
run_code3 = False
run_code4 = False

## USERS FUNCTIONS ##


#####################


## PROGRAM 1 ##

if run_code1 == True:
	print ""
	print "####################################"
	print "##### RUNNING 2D_calc_stack.py #####"
	print "####################################"
	print ""

	
#####################

## PROGRAM 2 ##

if run_code2 == True:
	print ""
	print "####################################"
	print "##### RUNNING 2D_calc_stack.py #####"
	print "#####	Using Galaxies	      #####"
	print "####################################"
	print ""
	
	use_gals = True
	code_num = 2

	## Individual Halos ##
	print ''
	print "## Working on Individual Halos ##"
	print ''
	gal_num = 100		# galaxies per halo for caustic surface
	run_num_ha = [0,100]	# run number for individual halos
	x_range,HaloID,M_crit200,R_crit200,R,V,INF_NFWMASS,DIA_CAUMASS,INF_CAU,INF_NFW,DIA_CAU = Calc.halo_gal_2D(gal_num,run_num_ha,code_num)

	## Ensemble Clusters ##
	print ''
	print '## Working on Ensemble Clusters ##'
	print ''
	scale_data = False	# For binning, do you want to scale r,v data or not?
	gal_num = 10
	bin_range = 10
	run_num_en = [0,10]	
	x_range,ENC_INF_NFWMASS,ENC_DIA_CAUMASS,ENC_INF_NFW,ENC_INF_CAU,ENC_DIA_CAU,ENC_R,ENC_V,ENC_M200,ENC_R200 = Calc.gal_stack_2D(gal_num,bin_range,run_num_en,code_num,scale_data)	

	ENC_INF_NFWMEAN = []
	ENC_DIA_CAUMEAN = []
	inf_nfw = copy(INF_NFW)
	dia_cau = copy(DIA_CAU)
	for j in range( run_num_en[1] ):
		inf_nfwmean = []
		dia_caumean = []
		for k in range(201):
			inf_nfwmean.append( mean( inf_nfw[j*bin_range:(j+1)*bin_range].T[k] ) )
			dia_caumean.append( mean( dia_cau[j*bin_range:(j+1)*bin_range].T[k] ) )
		ENC_INF_NFWMEAN.append( array(inf_nfwmean) )
		ENC_DIA_CAUMEAN.append( array(dia_caumean) )

	ENC_INF_NFWMEAN = array(ENC_INF_NFWMEAN)
	ENC_DIA_CAUMEAN =  array(ENC_DIA_CAUMEAN)


########################


## PROGRAM 4 ##
'''Single Halo Stacking, Self-Stacking'''
if run_code4 == True:
	print ""
	print "####################################"
	print "##### RUNNING 2D_calc_stack.py #####"
	print "#####	Using Galaxies	      #####"
	print "####################################"
	print ""

	## Constants ##
	line_num = 5	# Number of lines of sight

	## Flags ##


	## Initialization ##
	U = universal()
	
	## Program ##





















