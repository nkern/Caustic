# This program holds the classes for all 2D stacking

## IMPORT MODULES ##
from numpy import * 
import pyfits
from numpy.random import randint,uniform
from numpy.linalg import norm
from astStats import biweightScale
from flux_caustics_nideal import *
from scipy import weave
from scipy.weave import converters
from numpy import max as npmax
from scipy import *
import random

## CONSTANTS ##

h = 0.72	# Hubble constant/100.0

## INITIALIZATION ##

C = caustic()

######################

class universal():

	def load_halos(self,h,root):
		HaloID = loadtxt(''+str(root)+'/nkern/Documents/MDB_milliMil_halodata/Caustic/biglosclusters.csv', delimiter=',', dtype='string', usecols=(0,), unpack=True)
		HPX,HPY,HPZ,HVX,HVY,HVZ = loadtxt(''+str(root)+'/nkern/Documents/MDB_milliMil_halodata/100_halo_sample/biglosclusters.csv', delimiter=',', dtype='float', usecols=(9,10,11,12,13,14), unpack=True)
		SRAD,ESRAD,R_crit200,M_crit200,HVD,Z = loadtxt(''+str(root)+'/nkern/Documents/MDB_milliMil_halodata/Caustic/Millbig_concentrations.phys_phys.csv', delimiter=',', dtype='float', usecols=(1,2,5,7,9,12), unpack=True)	
		for l in xrange(len(HaloID)):	# Cosmological correction
			HPX[l],HPY[l],HPZ[l] = HPX[l]/(1+Z[l]),HPY[l]/(1+Z[l]),HPZ[l]/(1+Z[l])	
			# Fix weird SRAD values, if R_crit200/SRAD = Conc > 2, set SRAD=R_crit200/2
			if R_crit200[l]/h/SRAD[l] < 2.0:
				SRAD[l] = R_crit200[l]/h / 2.0

		#Ordering of halos in arrays is identical to biglosclusters' inherent ordering.
		return HaloID, R_crit200/h, M_crit200/h, HPX/h, HPY/h, HPZ/h, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z


	def sort_halos(self, HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z):
		sort = argsort(M_crit200)[::-1]	#Sorted by M_crit200 descending order
		HaloID = HaloID[sort]
		M_crit200 = M_crit200[sort]
		R_crit200 = R_crit200[sort]
		HPX = HPX[sort]
		HPY = HPY[sort]
		HPZ = HPZ[sort]
		HVX = HVX[sort]
		HVY = HVY[sort]
		HVZ = HVZ[sort]
		HVD = HVD[sort]
		SRAD = SRAD[sort]
		ESRAD = ESRAD[sort]
		Z = Z[sort]

		return HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z

	def bin_medcalc(self,bin_data):
		bin_value = median(copy(bin_data))
		return bin_value
	def bin_meancalc(self,bin_data):
		bin_value = mean(copy(bin_data))
		return bin_value

	def bin_mass_calc(self,M_crit200,j,bin_range):
		''' Calculates a particular Bin's Mass using a median'''
		m200 = copy(M_crit200)
		bin_mass = median(m200[j-bin_range+1:j+1])
		return bin_mass

	def enc_bin_table(self,M_crit200,R_crit200,SRAD,ESRAD,HVD,halo_num,bin_range):
		'''Bins table value data using a mean, constant w/ varying line of sight'''
		n = halo_num/bin_range
		enc_m200 = zeros(n)
		enc_r200 = zeros(n)
		enc_srad = zeros(n)
		enc_esrad = zeros(n)
		enc_hvd = zeros(n)
		for j in xrange(n):
			enc_m200[j] = median(M_crit200[j*bin_range:j*bin_range+bin_range])
			enc_r200[j] = mean(R_crit200[j*bin_range:j*bin_range+bin_range])
			enc_srad[j] = mean(SRAD[j*bin_range:j*bin_range+bin_range])
			enc_esrad[j] = mean(ESRAD[j*bin_range:j*bin_range+bin_range])
			enc_hvd[j] = mean(HVD[j*bin_range:j*bin_range+bin_range])
		return enc_m200,enc_r200,enc_srad,enc_esrad,enc_hvd

	def enc_value_medcalc(self,data,j,bin_range):
		''' Calculates a particular Bin's variable value using a median'''
		data2 = copy(data)
		enc_value = median(data2[j-bin_range+1:j+1])
		return enc_value


	def Pick_pos(self,halo_p):
	        '''Picks a random position for the observer a given distance away from the center'''
		x = random.uniform(-1,1)
		y = random.uniform(-1,1)
		z = random.uniform(-1,1)
		unit = array([x,y,z])/(x**2+y**2+z**2)**(.5)
		# move the position randomly 30Mpc away
	        return halo_p + 30*unit

	def line_of_sight(self,gal_p,gal_v,halo_p,halo_v,H0,c):
		'''Line of Sight Calculations'''
		# Pick Position
		new_pos = self.Pick_pos(halo_p) 

		# New Halo Information
		halo_dist = ((halo_p[0]-new_pos[0])**2 + (halo_p[1]-new_pos[1])**2 + (halo_p[2]-new_pos[2])**2)**0.5
		halo_pos_unit = array([halo_p[0]-new_pos[0],halo_p[1]-new_pos[1],halo_p[2]-new_pos[2]]) / halo_dist
		halo_vlos = dot(halo_pos_unit, halo_v)

		# New Galaxy Information
		gal_p = array(gal_p)
		gal_v = array(gal_v)
		gal_dist = ((gal_p[0]-new_pos[0])**2 + (gal_p[1]-new_pos[1])**2 + (gal_p[2]-new_pos[2])**2)**0.5
		gal_vlos = zeros(gal_dist.size)
		gal_pos_unit = zeros((3,gal_dist.size))	#vector from new_p to gal	
		n = gal_dist.size
		# Line of sight
		code = """
		int u,w;
		for (u=0;u<n;++u){
		for(w=0;w<3;++w){
		gal_pos_unit(w,u) = (gal_p(w,u)-new_pos(w))/gal_dist(u);
		}
		gal_vlos(u) = gal_pos_unit(0,u)*gal_v(0,u)+gal_pos_unit(1,u)*gal_v(1,u)+gal_pos_unit(2,u)*gal_v(2,u);
		}
		"""
		fast = weave.inline(code,['gal_pos_unit','n','gal_dist','gal_vlos','gal_v','new_pos','gal_p'],type_converters=converters.blitz,compiler='gcc')
		angles = arccos(dot(halo_pos_unit,gal_pos_unit))
		r = angles*halo_dist
		v_pec = gal_vlos-halo_vlos*dot(halo_pos_unit,gal_pos_unit)
		z_clus_cos = H0*halo_dist/c
		z_clus_pec = 0#halo_vlos/c
		z_clus_obs = (1+z_clus_pec)*(1+z_clus_cos)-1
		z_gal_cos = H0*gal_dist/c
		z_gal_pec = gal_vlos/c
		z_gal_obs = (1+z_gal_pec)*(1+z_gal_cos)-1
		v = c*(z_gal_obs-z_clus_obs)/(1+z_clus_obs)
		#gal_vdisp3d[i] = np.sqrt(astStats.biweightScale(gal_v[0][np.where(gal_radius<=HaloR200[i])]-Halo_V[0],9.0)**2+astStats.biweightScale(gal_v[1][np.where(gal_radius<=HaloR200[i])]-Halo_V[1],9.0)**2+astStats.biweightScale(gal_v[2][np.where(gal_radius<=HaloR200[i])]-Halo_V[2],9.0)**2)/np.sqrt(3)
		#print 'MY VELOCITY OF GALAXIES', gal_vdisp3d[i]
#		particle_vdisp3d[i] = HVD*np.sqrt(3)
#		gal_rmag_new = gal_abs_rmag# + 5*np.log10(gal_dist*1e6/10.0)

		return r, v 

	def shiftgapper(self,data,sort=True):
		npbin = 25
		gap_prev = 2000.0
		nbins = int(ceil(data[:,0].size/(npbin*1.0)))
		origsize = data[:,0].shape[0]
		data = data[argsort(data[:,0])] #sort by r to ready for binning
#		print 'NBINS FOR GAPPER = ', nbins
		for i in range(nbins):
#			print 'BEGINNING BIN:',str(i)
			databin = data[npbin*i:npbin*(i+1)]
#			print 'R BETWEEN', str(databin[:,0][0]),'and',str(databin[:,0][-1])
#			print 'DATA SIZE IN',databin[:,0].size
			datanew = None
			nsize = databin[:,0].size
			datasize = nsize-1
			if nsize > 5:
				while nsize - datasize > 0 and datasize >= 5:
#					print '    ITERATING'
					nsize = databin[:,0].size
					databinsort = databin[argsort(databin[:,1])] #sort by v
					f = (databinsort[:,1])[databinsort[:,1].size-int(ceil(databinsort[:,1].size/4.0))]-(databinsort[:,1])[int(ceil(databinsort[:,1].size/4.0))]
					gap = f/(1.349)
#					print '    GAP SIZE', str(gap)
					if gap < 500.0: break
					#if gap >= 2.0*gap_prev:
					#    gap = gap_prev
					#    print '   Altered gap = %.2f'%(gap)
					databelow = databinsort[databinsort[:,1]<=0]
					gapbelow =databelow[:,1][1:]-databelow[:,1][:-1]
					dataabove = databinsort[databinsort[:,1]>0]
					gapabove = dataabove[:,1][1:]-dataabove[:,1][:-1]
					try:
						if np.max(gapbelow) >= gap: vgapbelow = where(gapbelow >= gap)[0][-1]
						else: vgapbelow = -1
#						print 'MAX BELOW GAP',np.max(gapbelow)
						try: 
							datanew = append(datanew,databelow[vgapbelow+1:],axis=0)
						except:
							datanew = databelow[vgapbelow+1:]
					except ValueError:
						pass
					try:
						if np.max(gapabove) >= gap: vgapabove = where(gapabove >= gap)[0][0]
						else: vgapabove = 99999999
#						print 'MAX ABOVE GAP',np.max(gapabove)
						try: 
							datanew = append(datanew,dataabove[:vgapabove+1],axis=0)
						except:
							datanew = dataabove[:vgapabove+1]
					except ValueError:
						pass
					databin = datanew
					datasize = datanew[:,0].size
					datanew = None
#				print 'DATA SIZE OUT', databin[:,0].size
				if gap >= 500.0:
					gap_prev = gap
				else:
					gap_prev = 500.0
			try:
				datafinal = append(datafinal,databin,axis=0)
			except:
				datafinal = databin
#		print 'GALAXIES CUT =',str(origsize-datafinal[:,0].size)
		self.gal_cut = str(origsize-datafinal[:,0].size)
		if sort == True:	# re sort by mags
			datafinal = datafinal[argsort(datafinal[:,2])]	
		return datafinal



###### INITIALIZATION OF UNIVERSAL #######

U = universal()

class galaxies():

	def configure_galaxies(self,HaloID,h,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_limit,vlimit,R_crit200,HVD,halo_num,gal_num,root):
		MAGS = []
		Gal_V = []
		Gal_P = []
		Halo_P = []
		Halo_V = []
		gpapp = Gal_P.append
		gvapp = Gal_V.append

		ID = loadtxt(''+str(root)+'/nkern/Documents/MDB_milliMil_halodata/Caustic/cmiller.csv',delimiter=',',dtype='str',usecols=(0,),unpack=True)
		for haloid,k in zip(list(HaloID),list(xrange(halo_num))):
			mags, gpx3d, gpy3d, gpz3d, gvx3d, gvy3d, gvz3d = self.load_galaxies(ID,haloid,h,k,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_limit,R_crit200,HVD,root)
			mags = array(mags,float)
			halo_p = array([HPX[k],HPY[k],HPZ[k]]) 
			halo_v = array([HVX[k],HVY[k],HVZ[k]])
			gal_p = array([gpx3d,gpy3d,gpz3d],float)
			gal_v = array([gvx3d,gvy3d,gvz3d],float)
			
			MAGS.append(mags)
			Halo_P.append(halo_p)
			Halo_V.append(halo_v)			
			gpapp(gal_p)
			gvapp(gal_v)			

		Halo_P,Halo_V = array(Halo_P),array(Halo_V)
		MAGS = array(MAGS)		

		return Halo_P,Halo_V,Gal_P,Gal_V,MAGS


	def load_galaxies(self,ID,haloid,h,k,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_limit,R_crit200,HVD,root):
		
		f = pyfits.open(''+str(root)+'/giffordw/Millenium/30Mpchalos/'+str(haloid)+'.Guo30_2.fits')
		data = f[1].data
		z,gpx,gpy,gpz,gvx,gvy,gvz,mags = data.field(13),data.field(17),data.field(18),data.field(19),data.field(20),data.field(21),data.field(22),data.field(63)
		gpx,gpy,gpz = (gpx/(1+Z[k])/h),(gpy/(1+Z[k])/h),(gpz/(1+Z[k])/h)
		gvx,gvy,gvz = gvx-HVX[k],gvy-HVY[k],gvz-HVZ[k]
		mags = array(mags,float)
		BCG = where(gpx != HPX[k])	#Removing BCG from sample

		return mags[BCG], gpx[BCG], gpy[BCG], gpz[BCG], gvx[BCG], gvy[BCG], gvz[BCG]


	def scale_gals(self,r,v,r_crit200,hvd):
		r /= r_crit200
		v /= hvd
		return r, v

	def limit_gals(self,r,v,mags,r_crit200,gal_p,gal_v,r_limit,vlimit,gal_num,line_num,method_num,l):
		'''All data limiting is done here'''

		method = 'method'+str(method_num)+''	# Which method to stack?
		global samp_size
		samp_size = None
		# List out 3d data
		gpx3d,gpy3d,gpz3d = gal_p[0],gal_p[1],gal_p[2]
		gvx3d,gvy3d,gvz3d = gal_v[0],gal_v[1],gal_v[2]

		## Sort by Mag ##
		sort = argsort(mags)
		r,v,mags = r[sort],v[sort],mags[sort]
		gpx3d,gpy3d,gpz3d,gvx3d,gvy3d,gvz3d = gpx3d[sort],gpy3d[sort],gpz3d[sort],gvx3d[sort],gvy3d[sort],gvz3d[sort]

		## Limit Sample ##
		# Taking top brightest galaxies within r and v limits
		sample = where( (r<r_limit*r_crit200) & (v < vlimit) & (v > -vlimit) )
		r,v,mags = r[sample],v[sample],mags[sample]
		gpx3d,gpy3d,gpz3d,gvx3d,gvy3d,gvz3d = gpx3d[sample],gpy3d[sample],gpz3d[sample],gvx3d[sample],gvy3d[sample],gvz3d[sample]

		#####################	
		## Select Galaxies ##
		#####################
		if method == 'method0':
			'''Picking M Top Brightest, such that gal_num (or N) particles within r200'''
			### Build Ensemble (Put in a few extra to counter shiftgapper)
			within = where(r<r_crit200)[0]
			end = within[:gal_num][-1] + 1
			en_r,en_v,en_mags,en_gpx3d,en_gpy3d,en_gpz3d,en_gvx3d,en_gvy3d,en_gvz3d = r[0:end],v[0:end],mags[0:end],gpx3d[0:end],gpy3d[0:end],gpz3d[0:end],gvx3d[0:end],gvy3d[0:end],gvz3d[0:end] 
			### Build LOS (again, put a few extra in for shiftgapper)
			end = within[:gal_num + gal_num/5][-1] + 1
			rt,vt,magst,gpx3dt,gpy3dt,gpz3dt,gvx3dt,gvy3dt,gvz3dt = U.shiftgapper(vstack((r[:end],v[:end],mags[:end],gpx3d[:end],gpy3d[:end],gpz3d[:end],gvx3d[:end],gvy3d[:end],gvz3d[:end])).T).T
			within = where(rt<r_crit200)[0]
			end = within[:gal_num][-1] + 1
			r,v,mags,gpx3d,gpy3d,gpz3d,gvx3d,gvy3d,gvz3d = rt[:end],vt[:end],magst[:end],gpx3dt[:end],gpy3dt[:end],gpz3dt[:end],gvx3dt[:end],gvy3dt[:end],gvz3dt[:end]

		if method == 'method1':
			'''Randomly choosing M galaxies until gal_num galaxies are within r200'''
			sam_size = gal_num*25
			r,v,mags,gpx3d,gpy3d,gpz3d,gvx3d,gvy3d,gvz3d = r[:sam_size],v[:sam_size],mags[:sam_size],gpx3d[:sam_size],gpy3d[:sam_size],gpz3d[:sam_size],gvx3d[:sam_size],gvy3d[:sam_size],gvz3d[:sam_size]
			samp_size = len(r)
			while True:
				end = gal_num*1.5
				rando = randint(0,len(r),end)
				within = where(r[rando]<=r_crit200)[0]
				if len(within) < gal_num:
					end += 10	
				else:
					break
			### Build Ensemble
			rt,vt,magst,gpx3dt,gpy3dt,gpz3dt,gvx3dt,gvy3dt,gvz3dt = U.shiftgapper(vstack((r[rando],v[rando],mags[rando],gpx3d[rando],gpy3d[rando],gpz3d[rando],gvx3d[rando],gvy3d[rando],gvz3d[rando])).T).T	# decided to shiftgapper before stacking and los
			within = where(rt<r_crit200)[0]
			end = within[:gal_num+1][-1] + 1	# Leaving 1 extra gal per los for ensemble shiftgapper space
			en_r,en_v,en_mags,en_gpx3d,en_gpy3d,en_gpz3d,en_gvx3d,en_gvy3d,en_gvz3d = rt[0:end],vt[0:end],magst[0:end],gpx3dt[0:end],gpy3dt[0:end],gpz3dt[0:end],gvx3dt[0:end],gvy3dt[0:end],gvz3dt[0:end]
			### Build LOS
			end = within[:gal_num][-1] + 1
			r,v,mags,gpx3d,gpy3d,gpz3d,gvx3d,gvy3d,gvz3d = rt[:end],vt[:end],magst[:end],gpx3dt[:end],gpy3dt[:end],gpz3dt[:end],gvx3dt[:end],gvy3dt[:end],gvz3dt[:end]

		if method == 'method2':
			### Build Ensemble 
			within = where(r<r_crit200)[0]
			end = within[:gal_num*line_num + line_num][-1] + 1
			en_r,en_v,en_mags,en_gpx3d,en_gpy3d,en_gpz3d,en_gvx3d,en_gvy3d,en_gvz3d = r[:end][l::line_num],v[:end][l::line_num],mags[:end][l::line_num],gpx3d[:end][l::line_num],gpy3d[:end][l::line_num],gpz3d[:end][l::line_num],gvx3d[:end][l::line_num],gvy3d[:end][l::line_num],gvz3d[:end][l::line_num]	
			### Build LOS
			end = within[:gal_num*line_num + 5*line_num][-1] + 1
			rt,vt,magst,gpx3dt,gpy3dt,gpz3dt,gvx3dt,gvy3dt,gvz3dt = U.shiftgapper(vstack((r[:end][l::line_num],v[:end][l::line_num],mags[:end][l::line_num],gpx3d[:end][l::line_num],gpy3d[:end][l::line_num],gpz3d[:end][l::line_num],gvx3d[:end][l::line_num],gvy3d[:end][l::line_num],gvz3d[:end][l::line_num])).T).T
			within = where(rt<r_crit200)[0]
			end = within[:gal_num][-1] + 1
			r,v,mags,gpx3d,gpy3d,gpz3d,gvx3d,gvy3d,gvz3d = rt[:end],vt[:end],magst[:end],gpx3dt[:end],gpy3dt[:end],gpz3dt[:end],gvx3dt[:end],gvy3dt[:end],gvz3dt[:end] 


		return r,v,mags,gpx3d,gpy3d,gpz3d,gvx3d,gvy3d,gvz3d,en_r,en_v,en_mags,en_gpx3d,en_gpy3d,en_gpz3d,en_gvx3d,en_gvy3d,en_gvz3d 


	def bin_clusters(self,ENC_R,ENC_V,ENC_MAG,ENC_VDISP,ENC_R200,ENC_M200,ENC_SRAD,ENC_ESRAD,ENC_GPX3D,ENC_GPY3D,ENC_GPZ3D,ENC_GVX3D,ENC_GVY3D,ENC_GVZ3D,LINE_VDISP,Gal_P,Gal_V,Halo_P,Halo_V,M_crit200,R_crit200,SRAD,ESRAD,MAGS,k,r_limit,vlimit,gal_num,line_num,H0,q,c,LINE_DIAMASS,LINE_INFMASS,LINE_DIACAU,LINE_DISSECT,root,beta,scale_data):
	
		#update bin range: list of halo ids who belong in the kth ensemble
		bin_range = arange(k*line_num,k*line_num+line_num,1,int)

		#Binning arrays
		enc_r = []
		enc_v = []
		enc_mag = []
		enc_gpx3d = []
		enc_gpy3d = []
		enc_gpz3d = []
		enc_gvx3d = []
		enc_gvy3d = []
		enc_gvz3d = []

		#Line of sight arrays
		line_diamass = []
		line_infmass = []
		line_dia_cau = []
		line_inf_cau = []
		line_inf_nfw = []
		line_vdisp = []
		line_r = []
		line_v = []
		line_mag = []
		line_dissect = []

		#Loop over binned halos
		for l in bin_range:
						
			#Line of Sight Calculation
			r,v = U.line_of_sight(Gal_P[l],Gal_V[l],Halo_P[l],Halo_V[l],H0,c)
	
			#Limit Data
			r,v,mags,gal_vdisp,gpx3d,gpy3d,gpz3d,gvx3d,gvy3d,gvz3d = self.limit_gals(r,v,MAGS[l],R_crit200[l],Gal_P[l],Gal_V[l],r_limit,vlimit,gal_num,line_num,l)
			line_dissect.append( len(r) )

			# Append LOS RV arrays
			line_r.append(r)
			line_v.append(v)
			line_mag.append(mags)		
	
			#Scale Data
			if scale_data == True:
				r,v = self.scale_gals(r,v,R_crit200[l],gal_vdisp)

			# Do Mass Estimation for each Line of Sight 
			x_range,line_diamass,line_infmass,line_dia_cau,line_inf_cau,line_inf_nfw = self.kernel_caustic_masscalc(r,v,line_diamass,line_infmass,line_dia_cau,line_inf_cau,line_inf_nfw,R_crit200[l],M_crit200[l],SRAD[l],ESRAD[l],gal_vdisp,r_limit,vlimit,H0,q,k,root,beta,l=l)
			
			enc_r.extend(r)
			enc_v.extend(v)
			enc_mag.extend(mags)
			enc_gpx3d.extend(gpx3d)
			enc_gpy3d.extend(gpy3d)
			enc_gpz3d.extend(gpz3d)
			enc_gvx3d.extend(gvx3d)
			enc_gvy3d.extend(gvy3d)
			enc_gvz3d.extend(gvz3d)

			line_vdisp.append(gal_vdisp)

		# Shift Gapper Method to remove interlopers
		enc_r,enc_v,enc_mag,enc_gpx3d,enc_gpy3d,enc_gpz3d,enc_gvx3d,enc_gvy3d,enc_gvz3d = U.shiftgapper(vstack((enc_r,enc_v,enc_mag,enc_gpx3d,enc_gpy3d,enc_gpz3d,enc_gvx3d,enc_gvy3d,enc_gvz3d)).T).T

		# Calculated or Average Ensemble Properties
		enc_vdisp = biweightScale(copy(enc_v)[where( copy(enc_r)<R_crit200[l] )[0]],9.0)
		ENC_R200.append(U.bin_meancalc(R_crit200[bin_range]))
		ENC_M200.append(U.bin_medcalc(M_crit200[bin_range]))
		ENC_SRAD.append(U.bin_meancalc(SRAD[bin_range]))
		ENC_ESRAD.append(U.bin_meancalc(ESRAD[bin_range]))

		#Ensemble Arrays
		ENC_R.append(enc_r)
		ENC_V.append(enc_v)
		ENC_MAG.append(enc_mag)
		ENC_VDISP.append(enc_vdisp)
		ENC_GPX3D.append(enc_gpx3d)
		ENC_GPY3D.append(enc_gpy3d)
		ENC_GPZ3D.append(enc_gpz3d)
		ENC_GVX3D.append(enc_gvx3d)
		ENC_GVY3D.append(enc_gvy3d)
		ENC_GVZ3D.append(enc_gvz3d)

		#Line of Sight Arrays
		LINE_DIAMASS.append(line_diamass)
		LINE_INFMASS.append(line_infmass)
		LINE_DIACAU.append(line_dia_cau)
		LINE_VDISP.append(line_vdisp)
		LINE_DISSECT.append(line_dissect)

		return ENC_R,ENC_V,ENC_MAG,ENC_VDISP,ENC_R200,ENC_M200,ENC_SRAD,ENC_ESRAD,ENC_GPX3D,ENC_GPY3D,ENC_GPZ3D,ENC_GVX3D,ENC_GVY3D,ENC_GVZ3D,LINE_VDISP,LINE_DIAMASS,LINE_INFMASS,LINE_DIACAU,LINE_DISSECT


	def kernel_caustic_masscalc(self,enc_r,enc_v,ENC_DIA_MASS,ENC_INF_MASS,ENC_DIA_CAU,ENC_INF_CAU,ENC_INF_NFW,enc_r200,enc_m200,enc_srad,enc_esrad,enc_vdisp,r_limit,vlimit,H0,q,k,root,beta,l=None):
	
		enc_r,enc_v = array(enc_r),array(enc_v)
		r = append(enc_r,enc_r)
		v = append(enc_v,-1*enc_v)

		xbeta,abeta = loadtxt(''+str(root)+'/nkern/Documents/MDB_milliMil_halodata/Caustic/average_betaprofile.tab',dtype='float',usecols=(0,1),unpack=True)

		print ''
		print '## Working on Ensemble #'+str(k)+''
		if l != None:
			print '## Line of Sight #'+str(l)+''
		print '_______________________________________'
		fit = polyfit((xbeta*enc_r200)[xbeta<4],abeta[xbeta<4],6)			
		res_array = array([1])
		img_tot = 0
		img_grad_tot = 0
		img_inf_tot = 0
		for u in xrange(res_array.size):
			x_range,y_range,img,img_grad,img_inf = C.gaussian_kernel(r,v,enc_r200,normalization=H0,scale=q,res=200,adj=res_array[u],see=False)
			img_tot += img/npmax(img)
			img_grad_tot += img_grad/npmax(img_grad)
			img_inf_tot += img_inf/npmax(img_inf)

		## Define Beta ##
#		beta = fit[0]*x_range**6+fit[1]*x_range**5+fit[2]*x_range**4+fit[3]*x_range**3+fit[4]*x_range**2+fit[5]*x_range+fit[6]	
		beta = np.zeros(x_range.size) + beta
	
		## Caustic Surface Estimation ##
		maxv = enc_r200*H0*sqrt(200)+500

		## INFLECTION TECHNIQUE (MPHI) ##
		Anew,threesig,dens_norm,e_dens_norm,srad,e_srad = C.level_search2(r,v,r,v,r,x_range,y_range,img_tot,img_inf_tot,H0*q,enc_r200,r_limit,maxv,beta,enc_srad,enc_esrad,use_vdisp=enc_vdisp,bin=k+1)
		vdispersion = threesig/3.5

		## DIAFERIO TECHNIQUE (CAUSTIC) ##
		AnewD,threesigD,dens_normD,e_dens_normD,sradD,e_sradD = C.level_search(r,v,r,v,r,x_range,y_range,img_tot,H0*q,enc_r200,r_limit,maxv,beta,enc_srad,enc_esrad,use_vdisp=enc_vdisp,bin=k+1)

		## Mass Calculation ##
		# 1.) Using Diaferio Technique, pick out caustic, use caustic via massprofile to find mass, dia_caumass
		# 2.) Using Inflection technique find caustic, fit to an nfw, find mass via dens_norm, inf_nfwmass

		massprofile,integrand = C.masscalc(x_range,abs(C.Ar_final),enc_r200,enc_m200,vdispersion,beta=beta,conc=enc_r200/enc_srad)
		massprofileD,integrandD = C.masscalc(x_range,abs(C.Ar_finalD),enc_r200,enc_m200,vdispersion,beta=beta,conc=enc_r200/enc_srad)

		inf_nfwmass = 4*pi*dens_norm*(srad)**3*(log(1+enc_r200/srad)-enc_r200/srad/(1+enc_r200/srad))
		dia_caumass = massprofileD[where(x_range[x_range >= 0] < enc_r200)[0][-1]]

		
		ENC_DIA_MASS.append(dia_caumass)
		ENC_INF_MASS.append(inf_nfwmass)
		ENC_DIA_CAU.append(C.Ar_finalD)
		ENC_INF_CAU.append(C.Ar_final)
		ENC_INF_NFW.append(Anew)

		return x_range,ENC_DIA_MASS,ENC_INF_MASS,ENC_DIA_CAU,ENC_INF_CAU,ENC_INF_NFW


	################### SELF STACK ####################


	def self_stack_clusters(self,ENC_R,ENC_V,ENC_MAG,ENC_VDISP,ENC_GPX3D,ENC_GPY3D,ENC_GPZ3D,ENC_GVX3D,ENC_GVY3D,ENC_GVZ3D,LINE_VDISP,Gal_P,Gal_V,Halo_P,Halo_V,M_crit200,R_crit200,SRAD,ESRAD,MAGS,k,r_limit,vlimit,gal_num,line_num,method_num,H0,q,c,LINE_DIAMASS,LINE_INFMASS,LINE_DIACAU,LINE_R,LINE_V,LINE_MAG,root,beta):
		'''Not really binning data, but organizing self halo data to mimic binned data from before'''
		#Binning arrays
		enc_r = []
		enc_v = []
		enc_mag = []
		enc_gpx3d = []
		enc_gpy3d = []
		enc_gpz3d = []
		enc_gvx3d = []
		enc_gvy3d = []
		enc_gvz3d = []

		#Line of sight arrays
		line_diamass = []
		line_infmass = []
		line_dia_cau = []
		line_inf_cau = []
		line_inf_nfw = []
		line_vdisp = []
		line_r = []
		line_v = []
		line_mag = []

		for l in range(line_num):
						
			# Line of Sight Calculation
			r,v = U.line_of_sight(Gal_P[k],Gal_V[k],Halo_P[k],Halo_V[k],H0,c)
	
			# Limit Data
			r,v,mags,gpx3d,gpy3d,gpz3d,gvx3d,gvy3d,gvz3d,en_r,en_v,en_mags,en_gpx3d,en_gpy3d,en_gpz3d,en_gvx3d,en_gvy3d,en_gvz3d = self.limit_gals(r,v,MAGS[k],R_crit200[k],Gal_P[k],Gal_V[k],r_limit,vlimit,gal_num,line_num,method_num,l)

			# Build Ensemble Data (w/o gapping method per LOS)	
			enc_r.extend(en_r)
			enc_v.extend(en_v)
			enc_mag.extend(en_mags)
			enc_gpx3d.extend(en_gpx3d)
			enc_gpy3d.extend(en_gpy3d)
			enc_gpz3d.extend(en_gpz3d)
			enc_gvx3d.extend(en_gvx3d)
			enc_gvy3d.extend(en_gvy3d)
			enc_gvz3d.extend(en_gvz3d)

			# Calculate LOS HVD (after interloper removal)
			gal_count = len(where( r<R_crit200[k] )[0] )
			if gal_count <= 3:
				'''biweightScale freaks out w/ less than 3 data points'''
				gal_vdisp = std(copy(v)[where( r<R_crit200[k] )])
			if gal_count > 3:	# This is the best way to calculate vdisp
				gal_vdisp = biweightScale(copy(v)[where( r<R_crit200[k] )],9.0)

			# Running Mass Estimation on Line of Sight 
			x_range,line_diamass,line_infmass,line_dia_cau,line_inf_cau,line_inf_nfw = self.self_stack_kernel_caustic_masscalc(r,v,line_diamass,line_infmass,line_dia_cau,line_inf_cau,line_inf_nfw,R_crit200[k],M_crit200[k],SRAD[k],ESRAD[k],gal_vdisp,r_limit,vlimit,H0,q,k,root,beta,l,samp_size=samp_size)
	
			# Append LOS Data
			line_vdisp.append(gal_vdisp)
			line_r.append(r)
			line_v.append(v)
			line_mag.append(mags)

		# Shift Gapper Method to remove interlopers (for the ensemble)
		enc_r,enc_v,enc_mag,enc_gpx3d,enc_gpy3d,enc_gpz3d,enc_gvx3d,enc_gvy3d,enc_gvz3d = U.shiftgapper(vstack((enc_r,enc_v,enc_mag,enc_gpx3d,enc_gpy3d,enc_gpz3d,enc_gvx3d,enc_gvy3d,enc_gvz3d)).T).T
		# Reduce system to gal_num*line_num gals within r200
		within = where(enc_r<R_crit200[k])[0]
		end = within[:gal_num*line_num][-1] + 1
		enc_r,enc_v,enc_mag,enc_gpx3d,enc_gpy3d,enc_gpz3d,enc_gvx3d,enc_gvy3d,enc_gvz3d = enc_r[:end],enc_v[:end],enc_mag[:end],enc_gpx3d[:end],enc_gpy3d[:end],enc_gpz3d[:end],enc_gvx3d[:end],enc_gvy3d[:end],enc_gvz3d[:end] 
		# Calculate Ensemble HVD
		enc_vdisp = biweightScale(copy(enc_v)[where( enc_r<R_crit200[k] )],9.0)
		#Ensemble Arrays
		ENC_R.append(enc_r)
		ENC_V.append(enc_v)
		ENC_MAG.append(enc_mag)
		ENC_VDISP.append(enc_vdisp)
		ENC_GPX3D.append(enc_gpx3d)
		ENC_GPY3D.append(enc_gpy3d)
		ENC_GPZ3D.append(enc_gpz3d)
		ENC_GVX3D.append(enc_gvx3d)
		ENC_GVY3D.append(enc_gvy3d)
		ENC_GVZ3D.append(enc_gvz3d)

		#Line of Sight Arrays
		LINE_DIAMASS.append(line_diamass)
		LINE_INFMASS.append(line_infmass)
		LINE_DIACAU.append(line_dia_cau)
		LINE_VDISP.append(line_vdisp)
		LINE_R.append(line_r)
		LINE_V.append(line_v)
		LINE_MAG.append(line_mag)

		return ENC_R,ENC_V,ENC_MAG,ENC_VDISP,ENC_GPX3D,ENC_GPY3D,ENC_GPZ3D,ENC_GVX3D,ENC_GVY3D,ENC_GVZ3D,LINE_VDISP,LINE_DIAMASS,LINE_INFMASS,LINE_DIACAU,LINE_R,LINE_V,LINE_MAG


	def self_stack_kernel_caustic_masscalc(self,ENC_R,ENC_V,ENC_CAUMASS,ENC_INFMASS,ENC_CAUSURF,ENC_INFSURF,ENC_INFNFW,ENC_R200,ENC_M200,ENC_SRAD,ENC_ESRAD,ENC_VDISP,r_limit,vlimit,H0,q,k,root,beta,l=None,samp_size=None):
		
		ENC_R,ENC_V = array(ENC_R),array(ENC_V)
		R = append(ENC_R,ENC_R)
		V = append(ENC_V,-1*ENC_V)

		xbeta,abeta = loadtxt(''+str(root)+'/nkern/Documents/MDB_milliMil_halodata/Caustic/average_betaprofile.tab',dtype='float',usecols=(0,1),unpack=True)

		print ''
		print '## Working on Cluster #'+str(k)+''
		if l != None:
			print '## Line of Sight #'+str(l)+''
		print '_______________________________________'
		print 'galaxies within r200 =',where(ENC_R<ENC_R200)[0].size
		if samp_size != None:
			print 'sample size =',samp_size

		### Kernel Density Estimation ###
		fit = polyfit((xbeta*ENC_R200)[xbeta<4],abeta[xbeta<4],6)			
		res_array = array([1])
		img_tot = 0
		img_grad_tot = 0
		img_inf_tot = 0
		for u in xrange(res_array.size):
			x_range,y_range,img,img_grad,img_inf = C.gaussian_kernel(R,V,ENC_R200,normalization=H0,scale=q,res=200,adj=res_array[u],see=False)
			img_tot += img/npmax(img)
			img_grad_tot += img_grad/npmax(img_grad)
			img_inf_tot += img_inf/npmax(img_inf)

		### Define Beta ###
#			beta = fit[0]*x_range**6+fit[1]*x_range**5+fit[2]*x_range**4+fit[3]*x_range**3+fit[4]*x_range**2+fit[5]*x_range+fit[6]	
		beta = np.zeros(x_range.size) + beta
	
		### Caustic Surface Estimation ###
		maxv = ENC_R200*H0*sqrt(200)+500

		# Inflection Technique to find Surface
		Anew,threesig,dens_norm,e_dens_norm,srad,e_srad = C.level_search2(R,V,R,V,R,x_range,y_range,img_tot,img_inf_tot,H0*q,ENC_R200,r_limit,maxv,beta,ENC_SRAD,ENC_ESRAD,use_vdisp=ENC_VDISP,bin=k+1)
		vdispersion = threesig/3.5

		# Caustic Technique to find Surface
		AnewD,threesigD,dens_normD,e_dens_normD,sradD,e_sradD = C.level_search(R,V,R,V,R,x_range,y_range,img_tot,H0*q,ENC_R200,r_limit,maxv,beta,ENC_SRAD,ENC_ESRAD,use_vdisp=ENC_VDISP,bin=k+1)

		### Mass Calculation ###
		# Normal "Diaferio 1999 Mass" with fbeta = .65 (inherent setting is fbeta = .5)
		massprofile,integrand = C.masscalc(x_range,abs(C.Ar_finalD),ENC_R200,ENC_M200,vdispersion,beta=beta,conc=ENC_R200/ENC_SRAD,dstyle=True)
		cau_diamass = (.65/.5) * massprofile[where(x_range[x_range >= 0] < ENC_R200)[0][-1]]

		# Caustic Mass with Concentration parameter, no fbeta
		massprofile2,integrand2 = C.masscalc(x_range,abs(C.Ar_finalD),ENC_R200,ENC_M200,vdispersion,beta=beta,conc=ENC_R200/ENC_SRAD)
		cau_concmass = massprofile2[where(x_range[x_range >= 0] < ENC_R200)[0][-1]]

		# Inflection Mass using NFW
		inf_nfwmass = 4*pi*dens_norm*(srad)**3*(log(1+ENC_R200/srad)-ENC_R200/srad/(1+ENC_R200/srad))

		# Labeling Caustic Mass as Diaferio 1999 fbeta mass, in future, Caustic Mass is concentration mass
		ENC_CAUMASS.append(cau_diamass)	
		ENC_INFMASS.append(inf_nfwmass)
		ENC_CAUSURF.append(C.Ar_finalD)
		ENC_INFSURF.append(C.Ar_final)
		ENC_INFNFW.append(Anew)

		return x_range,ENC_CAUMASS,ENC_INFMASS,ENC_CAUSURF,ENC_INFSURF,ENC_INFNFW


class particles():

	def load_parts():
		return


###############		####################		##################
############### 	####################		##################


class stack_2D():

	def gal_stack_2D(self,gal_num,bin_range,run_num,code_num,scale_data):
		## This Program uses MDB galaxy data to stack halos and produce a mass estimate

		## DEFINE CONSTANTS ##

		h = 0.72 		# Hubble Constant / 100.0
		r_limit = 2		# Radius Limit of data in Mpc
		vlimit = 3500.0		# Velocity Disp limit in km/s
		H0 = h*100.0		# Hubble constant
		q = 10.0
		c = 300000.0
		halo_num = 100		# Total number of halos
		gal_num = gal_num	# Number of galaxies stacked per halo for en. clusters
		bin_range = bin_range	# Number of halos per ensemble cluster
		run_num = run_num	# Number of ensemble clusters to run caustic mass calc over

		## DEFINE FLAGS ##

		use_mems = False
		use_vdisp = True
		  
		## INITIALIZATION ##
		G = galaxies()
		C = caustic()
		U = universal()

		### PROGRAM ###

		print '...loading halos'
		HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.load_halos(h)
		HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.sort_halos(HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z)

		print '...loading galaxies'
		R, V, MAGS, GAL_VDISP, GPX, GPY = G.configure_galaxies(HaloID,h,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_limit,vlimit,R_crit200,HVD,halo_num,gal_num)

		if scale_data == True:
			print '...scaling galaxies'
			R, V = G.scale_gals(R,V,R_crit200,HVD)
		print '...binning ensembles'
		ENC_R, ENC_V, ENC_MAG, ENC_M200, ENC_R200, ENC_HVD, ENC_GAL_VDISP, ENC_SRAD, ENC_ESRAD = G.bin_data_mag(HaloID,R,V,MAGS,SRAD,ESRAD,M_crit200,R_crit200,HVD,GAL_VDISP,halo_num,bin_range,gal_num,scale_data)

		print '...caustic!'
		x_range,ENC_INF_NFWMASS,ENC_DIA_NFWMASS,ENC_INF_CAUMASS,ENC_DIA_CAUMASS,ENC_INF_MPROF,ENC_INF_NFW,ENC_INF_CAU,ENC_DIA_MPROF,ENC_DIA_NFW,ENC_DIA_CAU = G.kernel_caustic_masscalc(ENC_R,ENC_V,ENC_M200,ENC_R200,ENC_SRAD,ENC_ESRAD,ENC_HVD,ENC_GAL_VDISP,halo_num,bin_range,gal_num,H0,q,r_limit,run_num,use_vdisp,use_mems)

		return x_range,ENC_INF_NFWMASS,ENC_DIA_CAUMASS,ENC_INF_NFW,ENC_INF_CAU,ENC_DIA_CAU,ENC_R,ENC_V,ENC_M200,ENC_R200


	def part_stack_2D():

		return




	def halo_gal_2D(self,gal_num,run_num,code_num):
		''' Individual Halo, Projected Space Techniques'''

		## DEFINE FLAGS ##
		use_mems = False
		use_vdisp = True	# Feed table value hvd
		use_gals = True		# Use galaxies or particles?

		## DEFINE CONSTANTS ##
		h = 0.72 		# Hubble Constant / 100.0
		r_limit = 2		# Radius Limit of data in R_200
		H0 = h*100.0		# Hubble constant
		q = 10.0
		c = 300000.0
		bin_range = 1		# Technically it is working on ensemble code
		halo_num = 100		# Number of halos in sample
		run_num = run_num	# Number of halos to run program over, particles
		vlimit = 3500.0	
		gal_num = gal_num

		## INITIALIZATION ##
		C = caustic()
		G = galaxies()
		U = universal()

		## DEFINE FUNCTIONS ##

		#####	PROGRAM    ######

		print '...loading halos'
		HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.load_halos(h)
		HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.sort_halos(HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z)


		print '...loading galaxies'
		R, V, MAGS, GAL_VDISP, GPX, GPY = G.configure_galaxies(HaloID,h,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_limit,vlimit,R_crit200,HVD,halo_num,gal_num)


		print '...caustic!'
		x_range,INF_NFWMASS,DIA_NFWMASS,INF_CAUMASS,DIA_CAUMASS,INF_MPROF,INF_NFW,INF_CAU,DIA_MPROF,DIA_NFW,DIA_CAU = G.kernel_caustic_masscalc(R,V,M_crit200,R_crit200,SRAD,ESRAD,HVD,GAL_VDISP,halo_num,bin_range,gal_num,H0,q,r_limit,run_num,use_vdisp,use_mems)

		return x_range,HaloID,M_crit200,R_crit200,R,V,INF_NFWMASS,DIA_CAUMASS,INF_CAU,INF_NFW,DIA_CAU


	def halo_part_2D():
		
		return



