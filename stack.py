#	THIS CODE WORKS ON Umich/Astro/nkern account	

# Calling and gathering data on halos and their galaxies to be worked on with particle/galaxy stacking
# This uses Caustic Technique... perhaps outdated..

########## IMPORT MODULES ##########

from math import *
from numpy import *
from matplotlib.pyplot import *
from scipy import stats
from stack_class import *
from flux_caustics import *

########## FLAGS ##########

kde = True			#Do you want a kde plot?
use_gals = True			#Using galaxies, or particles?	(Don't use particles unless on flux, it takes too long)
use_big_sample = True		#Using millimillennium 100 halo sample or not? (else: 16 halo sample)
use_beta_profs	= False		#For 100 halo sample: do you want to use previously calculated beta profiles? (else: beta=constant) 

########## DEFINE CONSTANTS ##########

h = 0.72			# h = H0 / 100; determines cosmology
H0 = h*100.0			# H0
selection_range = 1.5		# coefficient of the r_crit200 value which determines data selection range
galnumber = 40			# number of galaxies (or particles) used for caustic surface estimation
Ob_distance = 100.0		# distance from halo for each line of sight calculation, in Mpc
q = 10.0			# parameter...		
c = 300000.0			# speed of light in km/s
lower_halo_range = 0		# lower bound for HaloID index in 100 sample mass bin		
upper_halo_range = 20 		# upper bound for HaloID index in 100 sample mass bin
beta_constant = 0.2		# if beta profiles aren't used, this value is used as a constant 

########## INITIALIZE CLASSES ##########

C = caustic()
MP = milli_particle()
MG = milli_galaxy()
BP = bigsample_particle()
BG = bigsample_galaxy()

########## LOAD AND CONFIGURE DATA ##########

print '..loading halos'

if use_big_sample == False:

	if use_gals == True:
		HaloID, m_crit200, r_crit200, halox, haloy, haloz, halovx, halovy, halovz, halo_p, halo_v, halo_vdisp, halo_conc = MG.load_halos(h)
	
	else:
		HaloID, m_crit200, r_crit200, halox, haloy, haloz, halovx, halovy, halovz, halo_vdisp, halo_p, halo_v  = MP.load_halos(h)
else:
	if use_gals == True:
		HaloID, m_crit200, r_crit200, halox, haloy, haloz, halovx, halovy, halovz, halo_p, halo_v, halo_vdisp, halo_conc, beta_profs = BG.load_halos(h,lower_halo_range,upper_halo_range,use_beta_profs)


print '..loading gals and peforming line of sight'

if use_big_sample == False:

	if use_gals == True:
		r,rvalues,vvalues,magvalues = MG.load_gals(h,HaloID,r_crit200,m_crit200,halox,haloy,haloz,halovx,halovy,halovz,halo_p,halo_v,selection_range,Ob_distance,halo_vdisp)		#rvalues and vvalues include total data sets and are standardized, r = radius from center of halo, not standardized.
	else:
		r,rvalues,vvalues = MP.load_particles(h,HaloID,r_crit200,m_crit200,halox,haloy,haloz,halovx,halovy,halovz,halo_p,halo_v,selection_range,Ob_distance,halo_vdisp)			

else:
	if use_gals == True:
		r, rvalues, vvalues, magvalues = BG.load_gals(h,HaloID,r_crit200,m_crit200,halox,haloy,haloz,halovx,halovy,halovz,halo_p,halo_v,selection_range,Ob_distance,halo_vdisp)		

print '..choosing data'

if use_big_sample == False:
	if use_gals == True:
		R, V, M = MG.choose_gals(rvalues,vvalues,magvalues,galnumber)
	else:
		R, V = MP.choose_particles(rvalues,vvalues,galnumber)

else:
	if use_gals == True:
		R, V, M = BG.choose_gals(rvalues,vvalues,magvalues,galnumber)


######## AVERAGE BACK TO VALUES WITH UNITS ##########

avg_r200 = sum(r_crit200)/len(r_crit200)
avg_vdisp = sum(halo_vdisp)/len(halo_vdisp)
avg_conc = sum(halo_conc)/len(halo_conc)
vlimit = 3500.0
rlimit = avg_r200*2
R = R*avg_r200
V = V*avg_vdisp


if use_gals == True:			
	Rvalues, Vvalues, Mvalues = MG.limit_sample(R,V,M,rlimit,vlimit,H0)
#else:
	#copy and paste from flux_los_part.py after line of sight calculation, vfix stuff...

########## KERNEL DENSITY ESTIMATION ##########

res_array = np.array([10,5,3,2,1])
img_tot = 0
for u in range(res_array.size):
	x_range,y_range,img = C.gaussian_kernel(R, V, r200=1.0, normalization=H0, scale=q, res=200, adj = res_array[u], see=False)
	img_tot += img/np.max(img)

if use_beta_profs == True:
	beta = BG.beta_calc(beta_profs,x_range)
else:
	beta = zeros(x_range.size) + beta_constant
	 


########## CAUSTIC SURFACE ESTIMATION ##########

maxv = avg_r200*H0*sqrt(200.0)+500.0
Anew, threesig = C.level_search(Rvalues,Vvalues,Rvalues,Vvalues,Mvalues,x_range,y_range,img_tot,H0*q,avg_r200,rlimit,maxv,beta,use_vdisp=avg_vdisp, use_mems=False)


########## CALCULATE MASS PROFILE ##########

massprofile = C.masscalc(x_range,abs(Anew),avg_r200,'realmass unknown',avg_vdisp,beta=beta, conc=avg_conc)

print massprofile[where(x_range[x_range>=0] < avg_r200)[0][-1]]

########## PLOTTING ##########

if use_gals == True:

	print '..plotting'
	
	fig = figure()
	ax1 = fig.add_subplot(211)
	ax1.plot(Rvalues,Vvalues,'k.', x_range, Anew,'b', [0,6], [0,0], 'r', markersize=3)
	xlabel('r / r200')
	ylabel('v / vdisp')

### KDE PLOT (OPTIONAL) ###

	if kde == True:
		Z = U.kde_plot(R,V)
		ax2 = fig.add_subplot(212)
		ax2.plot(R,V,'k.', markersize=3)	
		ax2.imshow(rot90(Z), cmap=cm.Spectral_r, extent=[U.rmin,U.rmax,U.vmin,U.vmax], aspect='auto')

	show()



