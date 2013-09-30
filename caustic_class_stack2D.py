## caustic_class_stack2D.py
'''
#########
-This program contains all the functions and classes used by caustic.mass_stack2D.py, which uses the caustic technique to estimate the mass of galaxy clusters, after applying a stacking method to create an ensemble.
-Parts of the program are either modelled or directly taken from D. Gifford's CausticMass.py code. 
-This is the most up-to-date stacking code.
#########

bin_stack:
	functions:
	attributes:

self_stack:
	functions:
	attributes:

caustic:
	functions:
	attributes:
'''

## Import Modules ##

class selfstack(object):

	def __init__(self,varib):
		''' Adding permanent program variables to class namespace'''
		self.__dict__.update(varib)


	def self_stack_clusters(self,HaloID,HaloData,Halo_P,Halo_V,Gal_P,Gal_V,Gal_Mags,k):
		''' Building Ensemble Cluster and Calculating Property Statistics '''
		# Unpack HaloData array 
		M_crit200,R_crit200,Z,SRAD,ESRAD,HVD = HaloData
		# Loop over lines of sight
		for l in range(self.line_num):
			if self.light_cone == True:
				# Configure RA, DEC and Z into cluster-centric radius and velocity
				pass
			else:
				# Line of sight calculation if necessary
				r, v = U.line_of_sight(Gal_P[k],Gal_V[k],Halo_P[k],Halo_V[k])

			# Limit Data
			r, v, mags = U.limit_gals(r,v,Gal_Mags[k],R_crit200[k])					



	


class binstack(object):
	def __init__(self,varib):
		''' Adding permanent program variables to class namespace'''
		self.__dict__.update(varib)


class caustic(object):
	def __init__(self,varib):
		''' Adding permanent program variables to class namespace'''
		self.__dict__.update(varib)



