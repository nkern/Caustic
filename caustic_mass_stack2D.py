## caustic.mass_stack2D.py
'''
########
-This program uses the Caustic Technique to estimate the mass of galaxy clusters, after having applied a stacking technique to create an ensemble cluster.
-This code works in tandem with caustic.class_stack2D.py
-Parts of the program are either modelled are taken directly from D. Gifford's CausticMass.py code.
-This is the most up-to-date stacking code.
-In general, uppercase data arrays refer to all halos, and lower case data arrays refer to a specific halo
######## 
'''

## IMPORT MODULES ##
print '...importing modules'
import numpy as np
import pyfits
from numpy.random import randint,uniform
from numpy.linalg import norm
import matplotlib.pyplot as mp
import cosmolopy.distance as cd
import astStats
import scipy.ndimage as ndi
from scipy.interpolate import interp1d
import random
import sys

from caustic_class_stack2D import *
from caustic_universal_stack2D import *

## FLAGS ##
self_stack	= True			# Run self-stack or bin-stack
scale_data	= False			# Scale data by r200 and vdisp if True
use_flux	= True			# Using Flux if True, using Sophie if False
write_data 	= True			# Write Data to Result directories if True
light_cone	= False			# Input RA|DEC projection data if True, if False inputting x,y,z 3D data
one_ens		= True			# Only solve for one ensemble cluster if true 
					# This is generally the case when using an HPC

## CONSTANTS ##
c 		= 2.99792e5		# speed of light in km/s
h		= 1.0			# Hubble Constant, unitless
H0		= h*100.0		# Hubble Constant, km s-1 Mpc-1
q		= 10.0			# Scale of Gaussian Kernel Density Estimator
beta		= 0.2			# Velocity Anisotropy Beta parameter, if constant profile
fbeta		= 0.65			# fbeta value, see 'Diaferio 1999'
r_limit 	= 1.25			# Radius Cut Scaled by R200
v_limit		= 3500.0		# Velocity Cut in km/s
data_set	= 'Guo30_2'		# Data set to draw semi analytic data from

halo_num	= 100			# Total number of halos loaded
ens_num		= int(sys.argv[1])	# Number of Ensembles to solve for
gal_num		= int(sys.argv[2])	# Number of galaxies taken per line of sight
line_num	= int(sys.argv[3])	# Number of lines of sight to stack over
method_num	= int(sys.argv[5])	# Self Stacking Method to Use

if use_flux == True: 
	root=str('/nfs/christoq_ls')	# Change directory scheme if using flux or sophie
else: 
	root=str('/n/Christoq1')

if self_stack == True:							# Change Write Directory Depending on Parameters
	write_loc = 'ss_m'+str(method_num)+'_run'+str(sys.argv[4])	# Self Stack data-write location
else:
	write_loc = 'bs_m'+str(method_num)+'_run'+str(sys.argv[4]) # Bin Stack data-write location

# Make dictionary for all variables that remain constant throughout duration of program!
varib = {'c':c,'h':h,'H0':H0,'q':q,'beta':beta,'fbeta':fbeta,'r_limit':r_limit,'v_limit':v_limit,'data_set':data_set,'halo_num':halo_num,'ens_num':ens_num,'gal_num':gal_num,'line_num':line_num,'method_num':method_num,'write_loc':write_loc,'root':root,'self_stack':self_stack,'scale_data':scale_data,'use_flux':use_flux,'write_data':write_data,'light_cone':light_cone,'one_ens':one_ens}

## INITIALIZATION ##
U = universal(varib)
BS = binstack(varib)
SS = selfstack(varib)
C  = caustic(varib)

###################
##### PROGRAM #####
###################
U.print_separation('Running caustic_mass_stack2D.py')

U.print_separation('...Loading Halos',type=2)
# Load Halo Data
HaloID,HaloData = U.load_halos()
# Sort Halos by A Priori Known Descending Mass (Mass Critical 200)
HaloID,HaloData = U.sort_halos(HaloID,HaloData)

# Load Galaxy Data
U.print_separation('...Loading Galaxies',type=2)
Halo_P,Halo_V,Gal_P,Gal_V,Gal_Mags,HaloData = U.configure_galaxies(HaloID,HaloData)

# Solve for Ensemble number: ens_num 
j = 0
for k in np.array([ens_num]):

	# Build Ensemble
	if self_stack:
#		SS.self_stack_clusters(HaloID,HaloData,Halo_P,Halo_V,Gal_P,Gal_V,Gal_Mags,k)
		pass
	else:
		BS.bin_stack_clusters()


	# Caustic Technique
	

	j += 1












