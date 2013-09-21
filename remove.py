'''A program to delete flux output data that I don't want, must be fed job id'''

import os
import sys
from numpy import arange


def rm_directory():
	root = 'ss_run'
	for k in arange(100):
		os.remove('stack_data/'+root+'/halo_'+str(k)+'_constants.tab')
		os.remove('stack_data/'+root+'/halo_'+str(k)+'_linenum.tab')
		os.remove('stack_data/'+root+'/halo_'+str(k)+'_losprofile.tab')
		os.remove('stack_data/'+root+'/halo_'+str(k)+'_profiles.tab')
		os.remove('stack_data/'+root+'/halo_'+str(k)+'_RVdata.tab')	
	os.remove('stack_data/'+root+'/program_constants.tab')

	return

def rm_flux():
	
	#Remove every file
	for id in arange(89,92):
		for k in range(100):
			os.remove('SELF-STACK.o102381'+str(id)+'-'+str(k)+'')

	return



