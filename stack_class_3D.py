# This program holds the classes for all 3D stacking

from numpy import *
import pyfits
from numpy.random import randint
from astStats import biweightScale
from flux_caustics_ideal import *
from numpy import max as npmax

### CONSTANTS ###

h = 0.72	# Hubble constant/100.0

### INITIALIZATION ###

C = caustic()

######

class universal():

	def load_halos(self,h):
		HaloID = loadtxt('/n/Christoq1/nkern/Documents/MDB_milliMil_halodata/Caustic/biglosclusters.csv', delimiter=',', dtype='string', usecols=(0,), unpack=True)
		HPX,HPY,HPZ,HVX,HVY,HVZ = loadtxt('/n/Christoq1/nkern/Documents/MDB_milliMil_halodata/100_halo_sample/biglosclusters.csv', delimiter=',', dtype='float', usecols=(9,10,11,12,13,14), unpack=True)
		SRAD,ESRAD,R_crit200,M_crit200,HVD,Z = loadtxt('/n/Christoq1/nkern/Documents/MDB_milliMil_halodata/Caustic/Millbig_concentrations.phys_phys.csv', delimiter=',', dtype='float', usecols=(1,2,5,7,9,12), unpack=True)	
		for l in range(len(HaloID)):
			HPX[l],HPY[l],HPZ[l] = HPX[l]/(1+Z[l]),HPY[l]/(1+Z[l]),HPZ[l]/(1+Z[l])	

		#Ordering of halos in arrays is identical to biglosclusters' inherent ordering.
		return HaloID, R_crit200/h, M_crit200/h, HPX/h, HPY/h, HPZ/h, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z


	def sort_halos(self, HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z):
		sort = argsort(M_crit200)[::-1]	#Sorted by M_crit200 descending order, can change at whim	
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


	def bin_mass_calc(self,M_crit200,j,bin_range):
		m200 = copy(M_crit200)
		bin_mass = median(m200[j-bin_range+1:j+1])
		return bin_mass

	def enc_value_calc(self,data,j,bin_range):
		data2 = copy(data)
		enc_value = mean(data2[j-bin_range+1:j+1])
		return enc_value


## INITIALIZE CLASS universal ###

U = universal()

######

class galaxies:

	def configure_galaxies(self,HaloID,h,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_limit,R_crit200,HVD,halo_num):
		R = []		# Each element of these lists will hold all galaxy data for a halo.
		V = []		# Therefore, each list will have 100 elements, each element being an array.
		MAGS = []
		GPX = []
		GPY = []
		GPZ = []

		ID = loadtxt('/n/Christoq1/nkern/Documents/MDB_milliMil_halodata/Caustic/cmiller.csv',delimiter=',',dtype='str',usecols=(0,),unpack=True)
		for haloid,k in zip(list(HaloID),list(range(halo_num))):
			r, v, mags, gpx, gpy, gpz = self.load_galaxies2(ID,haloid,h,k,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_limit,R_crit200,HVD)
			R.append(r)
			V.append(v)
			MAGS.append(mags)
			GPX.append(gpx)		
			GPY.append(gpy)
			GPZ.append(gpz)	

		R = array(R)
		V = array(V)
		MAGS = array(MAGS)
		GPX,GPY,GPZ = array(GPX),array(GPY),array(GPZ)

		return R, V, MAGS, GPX, GPY, GPZ


	
	def load_galaxies1(self,haloid,h,k,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_limit,R_crit200,HVD):
		gpx, gpy, gpz, gvx, gvy, gvz, mags, type = loadtxt('/n/Christoq1/nkern/Documents/MDB_milliMil_halodata/100_halo_sample/deeper_gal_sample/'+haloid+'_galaxies.csv', delimiter=',', dtype='float', usecols=(4,5,6,7,8,9,12,3), unpack=True)

		gpx,gpy,gpz = gpx/(1+Z[k]),gpy/(1+Z[k]),gpz/(1+Z[k])

		r = sqrt( (gpx/h-HPX[k])**2 + (gpy/h-HPY[k])**2 + (gpz/h-HPZ[k])**2 ) 
		v = sqrt( (gvx-HVX[k])**2 + (gvy-HVY[k])**2 + (gvz-HVZ[k])**2 )
		
		sort = argsort(mags)		# sorted by starting with brighest
		mags = array(mags[sort])
		type = array(type[sort])
		r = array(r[sort])/R_crit200[k]		#Standardizing to Halo's individual properties 
		v = array(v[sort])/HVD[k]

		## LIMIT DATA ##	
		cut = where((r<=r_limit) & (v<=5000/HVD[k]))
		mags = mags[cut]
		type = type[cut]
		r = r[cut]
		v = v[cut]

		return r, v, mags, type	

	def load_galaxies2(self,ID,haloid,h,k,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_limit,R_crit200,HVD):

		
		IDmatch = where(ID==haloid)[0][0]
		f = pyfits.open('/n/Christoq1/MILLENNIUM/particles/t_'+str(IDmatch)+'_cmiller_guo.fits')
		data = f[1].data
		z,gpx,gpy,gpz,gvx,gvy,gvz,mags = data.field(13),data.field(17),data.field(18),data.field(19),data.field(20),data.field(21),data.field(22),data.field(63)
		

		gpx,gpy,gpz = ( gpx/(1+Z[k])/h - HPX[k] ),( gpy/(1+Z[k])/h - HPY[k] ),( gpz/(1+Z[k])/h - HPZ[k] )
		gvx,gvy,gvz = gvx-HVX[k],gvy-HVY[k],gvz-HVZ[k]

		r = sqrt( (gpx)**2 + (gpy)**2 + (gpz)**2 ) 
		v = sqrt( (gvx)**2 + (gvy)**2 + (gvz)**2 )


		sort = argsort(mags)		# sorted by starting with brighest
		mags = array(mags[sort])
		r = array(r[sort])/R_crit200[k]		#Standardizing to Halo's individual properties 
		v = array(v[sort])/HVD[k]
		gpx,gpy,gpz,gvx,gvy,gvz = array(gpx[sort])/sqrt(R_crit200[k]),array(gpy[sort])/sqrt(R_crit200[k]),array(gpz[sort])/sqrt(R_crit200[k]),array(gvx[sort]),array(gvy[sort]),array(gvz[sort])
		## LIMIT DATA ##	
		cut = where((r<=r_limit) & (v<=5000.0/HVD[k]) & (v!=0)) #cutting out BCG
		r,v,mags = r[cut],v[cut],mags[cut]
		gpx,gpy,gpz,gvx,gvy,gvz = gpx[cut],gpy[cut],gpz[cut],gvx[cut],gvy[cut],gvz[cut]

	
		return r, v, mags, gpx, gpy, gpz


	def bin_data_mag(self,HaloID,R,V,MAGS,SRAD,ESRAD,M_crit200,R_crit200,HVD,halo_num,bin_range,gal_num,GPX,GPY,GPZ):

		# Galaxies chosen by brightest magnitude

		### Stacking Halos (Creating Ensemble clusters) ###
		ENC_R = []	#2 dimensional lists that hold ensemble cluster data (each bin is an ensemble cluster)
		ENC_V = []
		ENC_MAG = []
		ENC_GPX = []
		ENC_GPY = []	
		ENC_GPZ = []

		Rt = []		#Temporary lists to hold each bin's data
		Vt = []		#Lists are cleared after every bin
		Mt = []
		GPXt = []
		GPYt = []	
		GPZt = []
		#####

		### Calculating Averaged Ensemble Cluster values ###
		ENC_R200 = zeros(halo_num/bin_range)
		ENC_M200 = zeros(halo_num/bin_range)
		ENC_HVD = zeros(halo_num/bin_range)
		ENC_SRAD = zeros(halo_num/bin_range)
		ENC_ESRAD = zeros(halo_num/bin_range)
		#####

		for j in range(halo_num):
			Rt.extend(R[j][0:gal_num])		#Skip the BCG b/c it has no velocity
			Vt.extend(V[j][0:gal_num])
			Mt.extend(MAGS[j][0:gal_num])
			GPXt.extend(GPX[j][0:gal_num])
			GPYt.extend(GPY[j][0:gal_num])
			GPZt.extend(GPZ[j][0:gal_num])
			if j !=	0 and (j+1) % bin_range == 0:
				ENC_M200[(j+1)/bin_range-1] = U.bin_mass_calc(M_crit200,j,bin_range)
				ENC_R200[(j+1)/bin_range-1] = U.enc_value_calc(R_crit200,j,bin_range)
				ENC_HVD[(j+1)/bin_range-1] = U.enc_value_calc(HVD,j,bin_range)
				ENC_SRAD[(j+1)/bin_range-1] = U.enc_value_calc(SRAD,j,bin_range)
				ENC_ESRAD[(j+1)/bin_range-1] = U.enc_value_calc(ESRAD,j,bin_range)

				ENC_R.append(array(Rt))
				ENC_V.append(array(Vt))
				ENC_MAG.append(array(Mt))
				ENC_GPX.append(array(GPXt))
				ENC_GPY.append(array(GPYt))
				ENC_GPZ.append(array(GPZt))

				Rt = []
				Vt = []
				Mt = []
				GPXt = []
				GPYt = []
				GPZt = []

		ENC_R = array(ENC_R)
		ENC_V = array(ENC_V)
		ENC_MAG = array(ENC_MAG)
		ENC_GPX,ENC_GPY,ENC_GPZ = array(ENC_GPX),array(ENC_GPY),array(ENC_GPZ)

		for k in range(halo_num/bin_range):		#Putting real values on data
			ENC_R[k] *= ENC_R200[k]
			ENC_V[k] *= ENC_HVD[k]
			ENC_GPX[k] *= sqrt(ENC_R200[k])
			ENC_GPY[k] *= sqrt(ENC_R200[k])
			ENC_GPZ[k] *= sqrt(ENC_R200[k])
		return ENC_R, ENC_V, ENC_MAG, ENC_M200, ENC_R200, ENC_HVD, ENC_SRAD, ENC_ESRAD, ENC_GPX, ENC_GPY, ENC_GPZ
			
	def kernel_caustic_masscalc(self,ENC_R,ENC_V,ENC_M200,ENC_R200,ENC_SRAD,ENC_ESRAD,ENC_HVD,halo_num,bin_range,gal_num,H0,q,r_limit,run_num,use_mems):
		#Dummy arrays that Caustic functions call but have no use..
		potential = ones(halo_num/bin_range)
		ENC_GVD = ones(halo_num/bin_range)
	
		Rt=[]		# Doubling 3d data to mimic projected data, need for levelsearch2()
		Vt=[]
		for j in range(len(ENC_R)):
			Rt.append(append(ENC_R[j],ENC_R[j]))
			Vt.append(append(ENC_V[j],-ENC_V[j]))
		ENC_R = array(Rt)
		ENC_V = array(Vt)
		
		xbeta,abeta = loadtxt('/n/Christoq1/nkern/Documents/MDB_milliMil_halodata/Caustic/average_betaprofile.tab',dtype='float',usecols=(0,1),unpack=True)
		ENC_INF_MPROF = []
		ENC_INF_NFW = []
		ENC_INF_CAU = []

		ENC_DIA_MPROF = []
		ENC_DIA_NFW = []
		ENC_DIA_CAU = []
	
		ENC_DIA_CAUMASS = []	#Diaferio technique to produce Caustic, Caustic to find mass
		ENC_INF_CAUMASS = []	#Inflection technique to produce Casutic, Caustic to find mass	
		ENC_DIA_NFWMASS = []	#Diaferio technique to produce Caustic, NFW fit to find mass....?? check w/ dan on levelsearch() nfw fit, double parameter???	
		ENC_INF_NFWMASS = []	#Inflection technique to produce Caustic, NFW fit to find mass

		### Kernel Density Estimation ###
		for k in arange(run_num[0],run_num[1]):		#Loop over bins
			print ''
			print 'WORKING ON ENSEMBLE CLUSTER','#'+str(k)+''
			fit = polyfit((xbeta*ENC_R200[k])[xbeta<4],abeta[xbeta<4],6)	
			res_array = array([1])
			img_tot = 0
			img_grad_tot = 0
			img_inf_tot = 0
			for u in range(res_array.size):
				x_range,y_range,img,img_grad,img_inf = C.gaussian_kernel(ENC_R[k],ENC_V[k],ENC_R200[k],normalization=H0,scale=q,res=200,adj = res_array[u],see=False)
				img_tot += img/npmax(img)
				img_grad_tot += img_grad/npmax(img_grad)
				img_inf_tot += img_inf/npmax(img_inf)
		### Define Beta ###
			beta = fit[0]*x_range**6 + fit[1]*x_range**5 + fit[2]*x_range**4 + fit[3]*x_range**3 + fit[4]*x_range**2 + fit[5]*x_range + fit[6]

		### Caustic Surface Estimation ###
			maxv = ENC_R200[k]*H0*sqrt(200)+500
			
			## INFLECTION TECHNIQUE (M_PHI) ##
			Anew,threesig,dens_norm,e_dens_norm,srad,e_srad,Ar_final = C.level_search2(ENC_R[k],ENC_V[k],ENC_R[k],ENC_V[k],ENC_R[k],x_range,y_range,img_tot,img_inf_tot,H0*q,ENC_R200[k],r_limit,maxv,beta,potential,ENC_SRAD[k],ENC_ESRAD[k],use_vdisp=ENC_HVD[k],use_mems=use_mems,bin=k)
			vdispersion = threesig/3.5

			## DIAFERIO TECHNIQUE (CAUSTIC) ##
			Anewd,threesigd,dens_normd,e_dens_normd,sradd,e_sradd,Ar_finald = C.level_search(ENC_R[k],ENC_V[k],ENC_R[k],ENC_V[k],ENC_R[k],x_range,y_range,img_tot,H0*q,ENC_R200[k],r_limit,maxv,beta,potential,ENC_SRAD[k],ENC_ESRAD[k],use_vdisp=ENC_HVD[k],use_mems=use_mems,bin=k)
			vdispersiond = threesigd/3.5

		### Mass Calculation ###

			# Ways to calculate mass:
			#	1.) Using Inflection technique, pick out caustic, use caustic to calculate mass via massprofile
			#	2.) Using Diaferio technique, pick out casutic, use caustic to calculate mass via massprofiled
			#	3.) Using Inflection Causitc, fit an NFW, use NFW to calculate mass via dens_norm
			#	4.) Using Diaferio technique, fit an NFW, use NFW to calculate mass via dens_normd

			massprofile,integrand = C.masscalc(x_range,abs(Ar_final),ENC_R200[k],ENC_M200[k],vdispersion,beta,conc=ENC_R200[k]/ENC_SRAD[k])
			massprofiled,integrandd = C.masscalc(x_range,abs(Ar_finald),ENC_R200[k],ENC_M200[k],vdispersiond,beta,conc=ENC_R200[k]/ENC_SRAD[k])

			inf_nfwmass = 4*pi*dens_norm*(srad)**3*(np.log(1+ENC_R200[k]/srad)-ENC_R200[k]/srad/(1+ENC_R200[k]/srad))
			dia_nfwmass = 4*pi*dens_normd*(sradd)**3*(np.log(1+ENC_R200[k]/sradd)-ENC_R200[k]/sradd/(1+ENC_R200[k]/sradd))

			inf_caumass = massprofile[where(x_range[x_range >= 0] < ENC_R200[k])[0][-1]]			
			dia_caumass = massprofiled[where(x_range[x_range >= 0] < ENC_R200[k])[0][-1]]

		### Data Appending ###
			ENC_INF_NFWMASS.append(inf_nfwmass)
			ENC_INF_CAUMASS.append(inf_caumass)
			ENC_DIA_NFWMASS.append(dia_nfwmass)
			ENC_DIA_CAUMASS.append(dia_caumass)

			ENC_INF_MPROF.append(massprofile)
			ENC_INF_NFW.append(Anew)
			ENC_INF_CAU.append(Ar_final)

			ENC_DIA_MPROF.append(massprofiled)
			ENC_DIA_NFW.append(Anewd)
			ENC_DIA_CAU.append(Ar_finald)


		ENC_INF_MPROF,ENC_DIA_MPROF = array(ENC_INF_MPROF),array(ENC_DIA_MPROF)	
		ENC_INF_NFW,ENC_DIA_NFW = array(ENC_INF_NFW),array(ENC_DIA_NFW)
		ENC_INF_CAU,ENC_DIA_CAU = array(ENC_INF_CAU),array(ENC_DIA_CAU)
		

		return x_range,ENC_INF_NFWMASS,ENC_DIA_NFWMASS,ENC_INF_CAUMASS,ENC_DIA_CAUMASS,ENC_INF_MPROF,ENC_INF_NFW,ENC_INF_CAU,ENC_DIA_MPROF,ENC_DIA_NFW,ENC_DIA_CAU

class particles():


	def configure_particles(self,HaloID,h,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_limit,R_crit200,HVD,halo_num,gal_num,run_num):

		
		R = []		# Each element of these lists will hold all galaxy data for a halo.
		V = []		# Therefore, each list will have 100 elements, each element being an array.
		PPX = []
		PPY = []
		PPZ = []
		for haloid,k in zip(list(HaloID[run_num[0]:run_num[1]]),list(arange(run_num[0],run_num[1]))):
			r,v,ppx,ppy,ppz = self.load_particles(haloid,h,k,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_limit,R_crit200,HVD,gal_num)
			R.append(r)
			V.append(v)
			PPX.append(ppx)
			PPY.append(ppy)
			PPZ.append(ppz)

		R = array(R)
		V = array(V)
		PPX = array(ppx)
		PPY = array(ppy)
		PPZ = array(ppz)
		return R, V, PPX, PPY, PPZ


	def load_particles(self,haloid,h,k,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_limit,R_crit200,HVD,gal_num):


		id = loadtxt('/n/Christoq1/MILLENNIUM/particles/cmiller.csv', dtype='str', delimiter=',', usecols=(0,), unpack=True)
		id = delete(id,0) 
		index = where(id==str(haloid))

		p = pyfits.open('/n/Christoq1/MILLENNIUM/particles/t'+str(index[0][0])+'_cmiller.dat.fits')
		data = p[1].data
		ppx = data.field(1)/h/(1+Z[k])
		ppy = data.field(2)/h/(1+Z[k])
		ppz = data.field(3)/h/(1+Z[k])
		pvx = data.field(4)/sqrt(1+Z[k])
		pvy = data.field(5)/sqrt(1+Z[k])
		pvz = data.field(6)/sqrt(1+Z[k])

		pvx,pvy,pvz = pvx-HVX[k],pvy-HVY[k],pvz-HVZ[k]

		r = sqrt( (ppx**2) + (ppy**2) + (ppz**2) ) 
		v = sqrt( pvx**2 + pvy**2 + pvz**2 )

		r = array(r)/R_crit200[k]
		v = array(v)/HVD[k]		

		## LIMIT AND SELECT DATA ##
		cut = where((r<=r_limit) & (v<=5000.0/HVD[k]))
		r = r[cut]
		v = v[cut]
		ppx,ppy,ppz = ppx[cut],ppy[cut],ppz[cut]
		pick = randint(0,r.size,gal_num)	## RANDOM NUMBER PARTICLE SELECTION
		r,v,ppx,ppy,ppz = r[pick],v[pick],ppx[pick],ppy[pick],ppz[pick]
		print 'done loading halo',k
		return r, v, ppx, ppy, ppz

	def bin_data(self,HaloID,R,V,SRAD,ESRAD,M_crit200,R_crit200,HVD,halo_num,bin_range,run_num,gal_num):


		### Stacking Halos (Creating Ensemble clusters) ###
		ENC_R = []	#2 dimensional lists that hold ensemble cluster data 
		ENC_V = []
		Rt = []		#Temporary lists to hold each bin's data
		Vt = []		#Lists are cleared after every bin
		#####

		### Calculating Averaged Ensemble Cluster values ###
		ENC_R200 = zeros((run_num[1]-run_num[0])/bin_range)
		ENC_M200 = zeros((run_num[1]-run_num[0])/bin_range)
		ENC_HVD = zeros((run_num[1]-run_num[0])/bin_range)
		ENC_SRAD = zeros((run_num[1]-run_num[0])/bin_range)
		ENC_ESRAD = zeros((run_num[1]-run_num[0])/bin_range)
		#####

		for j in range(run_num[1]-run_num[0]):
			Rt.extend(R[j])
			Vt.extend(V[j])
			if j !=	0 and (j+1) % bin_range == 0:
				ENC_M200[(j+1)/bin_range-1] = U.bin_mass_calc(M_crit200,j,bin_range)
				ENC_R200[(j+1)/bin_range-1] = U.enc_value_calc(R_crit200,j,bin_range)
				ENC_HVD[(j+1)/bin_range-1] = U.enc_value_calc(HVD,j,bin_range)
				ENC_SRAD[(j+1)/bin_range-1] = U.enc_value_calc(SRAD,j,bin_range)
				ENC_ESRAD[(j+1)/bin_range-1] = U.enc_value_calc(ESRAD,j,bin_range)

				ENC_R.append(array(Rt))
				ENC_V.append(array(Vt))

				Rt = []
				Vt = []

		ENC_R = array(ENC_R)
		ENC_V = array(ENC_V)


		for k in range((run_num[1]-run_num[0])/bin_range):		#Putting real values on data
			ENC_R[k] *= ENC_R200[k]
			ENC_V[k] *= ENC_HVD[k]
		return ENC_R, ENC_V, ENC_M200, ENC_R200, ENC_HVD, ENC_SRAD, ENC_ESRAD


	def kernel_caustic_masscalc(self,ENC_R,ENC_V,ENC_M200,ENC_R200,ENC_SRAD,ENC_ESRAD,ENC_HVD,halo_num,bin_range,gal_num,H0,q,r_limit,run_num,use_mems):
		#Dummy arrays that Caustic functions call but have no use..
		potential = ones(halo_num/bin_range)
		ENC_GVD = ones(halo_num/bin_range)
	
		Rt=[]		# Doubling 3d data to mimic projected data, need for levelsearch2()
		Vt=[]
		for j in range(len(ENC_R)):
			Rt.append(append(ENC_R[j],ENC_R[j]))
			Vt.append(append(ENC_V[j],-ENC_V[j]))
		ENC_R = array(Rt)
		ENC_V = array(Vt)

		xbeta,abeta = loadtxt('/n/Christoq1/nkern/Documents/MDB_milliMil_halodata/Caustic/average_betaprofile.tab',dtype='float',usecols=(0,1),unpack=True)

		ENC_INF_MPROF = []
		ENC_INF_NFW = []
		ENC_INF_CAU = []

		ENC_DIA_MPROF = []
		ENC_DIA_NFW = []
		ENC_DIA_CAU = []
	
		ENC_DIA_CAUMASS = []	#Diaferio technique to produce Caustic, Caustic to find mass
		ENC_INF_CAUMASS = []	#Inflection technique to produce Casutic, Caustic to find mass	
		ENC_DIA_NFWMASS = []	#Diaferio technique to produce Caustic, NFW fit to find mass....?? check w/ dan on levelsearch() nfw fit, double parameter???	
		ENC_INF_NFWMASS = []	#Inflection technique to produce Caustic, NFW fit to find mass

		### Kernel Density Estimation ###
		for k in range((run_num[1]-run_num[0])/bin_range):
			''' k refers to the halo's position in ENC_* arrays, run_num/bin_range'''
			print ''
			print 'WORKING ON ENSEMBLE CLUSTER','#'+str(k)+''
			fit = polyfit((xbeta*ENC_R200[k])[xbeta<4],abeta[xbeta<4],6)	
			res_array = array([1])
			img_tot = 0
			img_grad_tot = 0
			img_inf_tot = 0
			for u in range(res_array.size):
				x_range,y_range,img,img_grad,img_inf = C.gaussian_kernel(ENC_R[k],ENC_V[k],ENC_R200[k],normalization=H0,scale=q,res=200,adj = res_array[u],see=False)
				img_tot += img/npmax(img)
				img_grad_tot += img_grad/npmax(img_grad)
				img_inf_tot += img_inf/npmax(img_inf)
		### Define Beta ###
			beta = fit[0]*x_range**6 + fit[1]*x_range**5 + fit[2]*x_range**4 + fit[3]*x_range**3 + fit[4]*x_range**2 + fit[5]*x_range + fit[6]

		### Caustic Surface Estimation ###
			maxv = ENC_R200[k]*H0*sqrt(200)+500
			
			## INFLECTION TECHNIQUE (M_PHI) ##
			Anew,threesig,dens_norm,e_dens_norm,srad,e_srad,Ar_final = C.level_search2(ENC_R[k],ENC_V[k],ENC_R[k],ENC_V[k],ENC_R[k],x_range,y_range,img_tot,img_inf_tot,H0*q,ENC_R200[k],r_limit,maxv,beta,potential,ENC_SRAD[k],ENC_ESRAD[k],use_vdisp=ENC_HVD[k],use_mems=use_mems,bin=k)
			vdispersion = threesig/3.5

			## DIAFERIO TECHNIQUE (CAUSTIC) ##
			Anewd,threesigd,dens_normd,e_dens_normd,sradd,e_sradd,Ar_finald = C.level_search(ENC_R[k],ENC_V[k],ENC_R[k],ENC_V[k],ENC_R[k],x_range,y_range,img_tot,H0*q,ENC_R200[k],r_limit,maxv,beta,potential,ENC_SRAD[k],ENC_ESRAD[k],use_vdisp=ENC_HVD[k],use_mems=use_mems,bin=k)
			vdispersiond = threesigd/3.5

		### Mass Calculation ###

			# Ways to calculate mass:
			#	1.) Using Inflection technique, pick out caustic, use caustic to calculate mass via massprofile
			#	2.) Using Diaferio technique, pick out casutic, use caustic to calculate mass via massprofiled
			#	3.) Using Inflection Causitc, fit an NFW, use NFW to calculate mass via dens_norm
			#	4.) Using Diaferio technique, fit an NFW, use NFW to calculate mass via dens_normd

			massprofile,integrand = C.masscalc(x_range,abs(Ar_final),ENC_R200[k],ENC_M200[k],vdispersion,beta,conc=ENC_R200[k]/ENC_SRAD[k])
			massprofiled,integrandd = C.masscalc(x_range,abs(Ar_finald),ENC_R200[k],ENC_M200[k],vdispersiond,beta,conc=ENC_R200[k]/ENC_SRAD[k])

			inf_nfwmass = 4*pi*dens_norm*(srad)**3*(np.log(1+ENC_R200[k]/srad)-ENC_R200[k]/srad/(1+ENC_R200[k]/srad))
			dia_nfwmass = 4*pi*dens_normd*(sradd)**3*(np.log(1+ENC_R200[k]/sradd)-ENC_R200[k]/sradd/(1+ENC_R200[k]/sradd))

			inf_caumass = massprofile[where(x_range[x_range >= 0] < ENC_R200[k])[0][-1]]			
			dia_caumass = massprofiled[where(x_range[x_range >= 0] < ENC_R200[k])[0][-1]]

		### Data Appending ###
			ENC_INF_NFWMASS.append(inf_nfwmass)
			ENC_INF_CAUMASS.append(inf_caumass)
			ENC_DIA_NFWMASS.append(dia_nfwmass)
			ENC_DIA_CAUMASS.append(dia_caumass)

			ENC_INF_MPROF.append(massprofile)
			ENC_INF_NFW.append(Anew)
			ENC_INF_CAU.append(Ar_final)

			ENC_DIA_MPROF.append(massprofiled)
			ENC_DIA_NFW.append(Anewd)
			ENC_DIA_CAU.append(Ar_finald)


		ENC_INF_MPROF,ENC_DIA_MPROF = array(ENC_INF_MPROF),array(ENC_DIA_MPROF)	
		ENC_INF_NFW,ENC_DIA_NFW = array(ENC_INF_NFW),array(ENC_DIA_NFW)
		ENC_INF_CAU,ENC_DIA_CAU = array(ENC_INF_CAU),array(ENC_DIA_CAU)
		

		return x_range,ENC_INF_NFWMASS,ENC_DIA_NFWMASS,ENC_INF_CAUMASS,ENC_DIA_CAUMASS,ENC_INF_MPROF,ENC_INF_NFW,ENC_INF_CAU,ENC_DIA_MPROF,ENC_DIA_NFW,ENC_DIA_CAU



###############		####################		##################
############### 	####################		##################



class stack_3D():	# This class contains versions of the '3D_*.py' programs in function form

	def gal_stack_3D(self,bin_range,gal_num,run_num,code_num):	

		'''This function is the program 3D_gal_stack.py'''

		# last update: 1/24/13

		# This program takes 3D galaxy data of 100 halo sample from Gerard Lemson from the MDB
		# and stacks the data by mass bin and uses the M_Phi technique. Note:(each bin is an ensemble cluster)

		import cosmolopy.distance as cd

		## DEFINE CONSTANTS ##

		h = 0.72 		# Hubble Constant / 100.0
		r_limit = 2		# Radius Limit of data in terms of R_crit200
		H0 = h*100.0		# Hubble constant
		q = 10.0
		c = 300000.0
		cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'h':H0/100.0}
		cosmo = cd.set_omega_k_0(cosmo)
		halo_num = 100		# Total number of halos
		gal_num = gal_num	# Number of galaxies stacked per halo for en. clusters
		bin_range = bin_range	# Number of halos per ensemble cluster
		run_num = run_num	# Number of ensemble halos to run caustic mass est. over

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

		R, V, MAGS, GPX, GPY, GPZ = G.configure_galaxies(HaloID,h,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_limit,R_crit200,HVD,halo_num)

		print '...binning data'	
		# All variables beginning with 'ENC_' stand for ensemble cluster variable, similar to *bin variables

		ENC_R,ENC_V,ENC_MAG,ENC_M200,ENC_R200,ENC_HVD,ENC_SRAD,ENC_ESRAD,ENC_GPX,ENC_GPY,ENC_GPZ = G.bin_data_mag(HaloID,R,V,MAGS,SRAD,ESRAD,M_crit200,R_crit200,HVD,halo_num,bin_range,gal_num,GPX,GPY,GPZ)

		print '...running caustic'

		x_range,ENC_INF_NFWMASS,ENC_DIA_NFWMASS,ENC_INF_CAUMASS,ENC_DIA_CAUMASS,ENC_INF_MPROF,ENC_INF_NFW,ENC_INF_CAU,ENC_DIA_MPROF,ENC_DIA_NFW,ENC_DIA_CAU = G.kernel_caustic_masscalc(ENC_R,ENC_V,ENC_M200,ENC_R200,ENC_SRAD,ENC_ESRAD,ENC_HVD,halo_num,bin_range,gal_num,H0,q,r_limit,run_num,use_mems)

		print ''
		bias1 = mean( (ENC_M200[run_num[0]:run_num[1]]-ENC_INF_NFWMASS) / ENC_M200[run_num[0]:run_num[1]] )
		bias2 = mean( abs(ENC_M200[run_num[0]:run_num[1]]-ENC_INF_NFWMASS) / ENC_M200[run_num[0]:run_num[1]] )
		bias3 = mean( log(ENC_INF_NFWMASS/ENC_M200[run_num[0]:run_num[1]]) )	
		
		if code_num == 1:	
			return bias1, bias2, -1*bias3
		if code_num == 2:
			return x_range,ENC_INF_NFWMASS,ENC_DIA_NFWMASS,ENC_INF_CAUMASS,ENC_DIA_CAUMASS,ENC_INF_MPROF,ENC_INF_NFW,ENC_INF_CAU,ENC_DIA_MPROF,ENC_DIA_NFW,ENC_DIA_CAU,ENC_M200,ENC_R200,ENC_R,ENC_V,ENC_MAG



	def part_stack_3D(self,bin_range,gal_num,run_num,code_num):		

		''' This function is the program 3D_part_stack.py'''
		# This program takes 3D particle data of 100 halo sample from Gerard Lemson from the MDB
		# and stacks the data by mass bin and uses the M_Phi technique. Note:(each bin is an ensemble cluster)

		# last update: 1/29/13

		##########

		import cosmolopy.distance as cd

		## DEFINE CONSTANTS ##

		h = 0.72 		# Hubble Constant / 100.0
		r_limit = 2		# Radius Limit of data in terms of R_crit200
		H0 = h*100.0		# Hubble constant
		q = 10.0
		c = 300000.0
		cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'h':H0/100.0}
		cosmo = cd.set_omega_k_0(cosmo)
		halo_num = 100		# Total number of halos

		## DEFINE FLAGS ##

		use_mems = False
		use_vdisp = True

		## INITIALIZATION ##

		G = galaxies()
		P = particles()
		C = caustic()
		U = universal()

		### PROGRAM ###

		print '...loading halos'

		HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.load_halos(h)

		HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.sort_halos(HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z)

		print '...loading particles'

		R, V, PPX, PPY, PPZ = P.configure_particles(HaloID,h,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,r_limit,R_crit200,HVD,halo_num,gal_num,run_num)

		print '...binning data'
		# All variables beginning with 'ENC_' stand for ensemble cluster, same for *bin variables

		ENC_R,ENC_V,ENC_M200,ENC_R200,ENC_HVD,ENC_SRAD,ENC_ESRAD = P.bin_data(HaloID,R,V,SRAD,ESRAD,M_crit200,R_crit200,HVD,halo_num,bin_range,run_num,gal_num)

		print '...running caustic'

		x_range,ENC_INF_NFWMASS,ENC_DIA_NFWMASS,ENC_INF_CAUMASS,ENC_DIA_CAUMASS,ENC_INF_MPROF,ENC_INF_NFW,ENC_INF_CAU,ENC_DIA_MPROF,ENC_DIA_NFW,ENC_DIA_CAU = P.kernel_caustic_masscalc(ENC_R,ENC_V,ENC_M200,ENC_R200,ENC_SRAD,ENC_ESRAD,ENC_HVD,halo_num,bin_range,gal_num,H0,q,r_limit,run_num,use_mems)

		return x_range,ENC_INF_NFWMASS,ENC_DIA_NFWMASS,ENC_INF_CAUMASS,ENC_DIA_CAUMASS,ENC_INF_MPROF,ENC_INF_NFW,ENC_INF_CAU,ENC_DIA_MPROF,ENC_DIA_NFW,ENC_DIA_CAU,ENC_R,ENC_V,ENC_M200,ENC_R200 


	def halo_gal_3D(self,gal_num,run_num,code_num):		

		'''This function is the program 3D_halos.py'''

		# 3d galaxy m_phi code for non-stacked halos.

		# last update: 1/29/13
		
		############
		import cosmolopy.distance as cd

		## DEFINE FLAGS ##

		use_mems = False
		use_vdisp = True

		## DEFINE CONSTANTS ##

		h = 0.72 		# Hubble Constant / 100.0
		r_limit = 2		# Radius Limit of data in terms of R_crit200
		H0 = h*100.0		# Hubble constant
		q = 10.0
		c = 300000.0
		cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'h':H0/100.0}
		cosmo = cd.set_omega_k_0(cosmo)
		halo_num = 100		# Number of halos in sample
		run_num = run_num	# Number of halos to run program over, particles
		bin_range = 1		# Needed b/c technically it is working on ensemble code

		## DEFINE FUNCTIONS ##

		def load_gals(h,gal_num,HaloID,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,R_crit200):
			R = []
			V = []
			MAGS = []
			GPX = []
			GPY = []
			GPZ = []	
			ID = loadtxt('/n/Christoq1/nkern/Documents/MDB_milliMil_halodata/Caustic/cmiller.csv',delimiter=',',dtype='str',usecols=(0,),unpack=True)
			for haloid,k in zip(list(HaloID),list(range(halo_num))):
				IDmatch = where(ID==str(haloid))[0][0]
				f = pyfits.open('/n/Christoq1/MILLENNIUM/particles/t_'+str(IDmatch)+'_cmiller_guo.fits')
				data = f[1].data
				z,gpx,gpy,gpz,gvx,gvy,gvz,mags = data.field(13),data.field(17),data.field(18),data.field(19),data.field(20),data.field(21),data.field(22),data.field(63)
				
				gpx,gpy,gpz = ( gpx/(1+Z[k])/h - HPX[k] ),( gpy/(1+Z[k])/h - HPY[k] ),( gpz/(1+Z[k])/h - HPZ[k] )
				gvx,gvy,gvz = gvx-HVX[k],gvy-HVY[k],gvz-HVZ[k]
			
				r = sqrt( (gpx)**2 + (gpy)**2 + (gpz)**2 ) 
				v = sqrt( (gvx)**2 + (gvy)**2 + (gvz)**2 )

				sort = argsort(mags)		# sorted by descending magnitude
				mags = array(mags[sort])
				r = array(r[sort])		 
				v = array(v[sort])
				gpx,gpy,gpz,gvx,gvy,gvz = array(gpx[sort]),array(gpy[sort]),array(gpz[sort]),array(gvx[sort]),array(gvy[sort]),array(gvz[sort])
				## LIMIT DATA ##
				cut = where((r<=r_limit*R_crit200[k]) & (v<=5000.0) & (v!=0) )[0][0:gal_num]	#skip BCG, no V 
				
				r,v,mags = r[cut],v[cut],mags[cut]
				gpx,gpy,gpz,gvx,gvy,gvz = gpx[cut],gpy[cut],gpz[cut],gvx[cut],gvy[cut],gvz[cut]	

				R.append(r)
				V.append(v)
				MAGS.append(mags)
				GPX.append(gpx)
				GPY.append(gpy)
				GPZ.append(gpz)
				
			R = array(R)
			V = array(V)
			MAGS = array(MAGS)
			GPX = array(GPX)
			GPY = array(GPY)
			GPZ = array(GPZ)
			
			return R,V,MAGS,GPX,GPY,GPZ	

			
		## INITIALIZATION ##

		U = universal()
		P = particles()
		G = galaxies()
		C = caustic()

		### PROGRAM ###

		print '...loading halos'

		HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.load_halos(h)

		HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.sort_halos(HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z)

		print '...loading gals'
		R, V, MAGS, GPX, GPY, GPZ = load_gals(h,gal_num,HaloID,HPX,HPY,HPZ,HVX,HVY,HVZ,Z,R_crit200)

		print '...caustic!'
		x_range,INF_NFWMASS,DIA_NFWMASS,INF_CAUMASS,DIA_CAUMASS,INF_MPROF,INF_NFW,INF_CAU,DIA_MPROF,DIA_NFW,DIA_CAU = G.kernel_caustic_masscalc(R,V,M_crit200,R_crit200,SRAD,ESRAD,HVD,halo_num,bin_range,gal_num,H0,q,r_limit,run_num,use_mems)

		return x_range,INF_NFWMASS,DIA_NFWMASS,INF_CAUMASS,DIA_CAUMASS,INF_MPROF,INF_NFW,INF_CAU,DIA_MPROF,DIA_NFW,DIA_CAU,R,V,MAGS,M_crit200,R_crit200


	def halo_part_3D(self,gal_num,run_num,code_num):

		# 3d particle m_phi code for non-stacked halos.

		# last update: 1/29/13
		
		############
		import cosmolopy.distance as cd
		from numpy.random import randint
		## DEFINE FLAGS ##

		use_mems = False
		use_vdisp = True

		## DEFINE CONSTANTS ##

		h = 0.72 		# Hubble Constant / 100.0
		r_limit = 2		# Radius Limit of data in terms of R_crit200
		H0 = h*100.0		# Hubble constant
		q = 10.0
		c = 300000.0
		cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'h':H0/100.0}
		cosmo = cd.set_omega_k_0(cosmo)
		halo_num = 100		# Number of halos in sample
		run_num = run_num	# Number of halos to run program over, particles
		bin_range = 1		# Needed b/c technically working on ensemble code

		## DEFINE FUNCTIONS ##

		def load_parts(h,gal_num,r_limit,HaloID,HVX,HVY,HVZ,Z,R_crit200,HVD):
			R = []
			V = []
			PPX = []
			PPY = []
			PPZ = []

			for haloid,k in zip(list(HaloID[run_num[0]:run_num[1]]),list(arange(run_num[0],run_num[1]))):
				id = loadtxt('/n/Christoq1/MILLENNIUM/particles/cmiller.csv', dtype='str', delimiter=',', usecols=(0,), unpack=True)
				id = delete(id,0) 
				index = where(id==str(haloid))

				p = pyfits.open('/n/Christoq1/MILLENNIUM/particles/t'+str(index[0][0])+'_cmiller.dat.fits')
				data = p[1].data
				ppx = data.field(1)/h/(1+Z[k])
				ppy = data.field(2)/h/(1+Z[k])
				ppz = data.field(3)/h/(1+Z[k])
				pvx = data.field(4)/sqrt(1+Z[k])
				pvy = data.field(5)/sqrt(1+Z[k])
				pvz = data.field(6)/sqrt(1+Z[k])

				pvx,pvy,pvz = pvx-HVX[k],pvy-HVY[k],pvz-HVZ[k]

				r = sqrt( (ppx**2) + (ppy**2) + (ppz**2) ) 
				v = sqrt( pvx**2 + pvy**2 + pvz**2 )

				r = array(r)
				v = array(v)		

				## LIMIT AND SELECT DATA ##
				cut = where((r<=r_limit*R_crit200[k]) & (v<=5000.0))
				r = r[cut]
				v = v[cut]
				ppx,ppy,ppz = ppx[cut],ppy[cut],ppz[cut]
				pick = randint(0,r.size,gal_num)	## RANDOM NUMBER PARTICLE SELECTION
				r,v,ppx,ppy,ppz = r[pick],v[pick],ppx[pick],ppy[pick],ppz[pick]

				R.append(r)
				V.append(v)
				PPX.append(ppx)
				PPY.append(ppy)
				PPZ.append(ppz)				

				print 'done loading halo',k

			R = array(R)
			V = array(V)
			PPX = array(PPX)
			PPY = array(PPY)
			PPZ = array(PPZ)

			return R, V, PPX, PPY, PPZ	
			
			
		## INITIALIZATION ##

		U = universal()
		P = particles()
		G = galaxies()
		C = caustic()

		### PROGRAM ###

		print '...loading halos'

		HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.load_halos(h)

		HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z = U.sort_halos(HaloID, R_crit200, M_crit200, HPX, HPY, HPZ, HVX, HVY, HVZ, HVD, SRAD, ESRAD, Z)

		print '...loading particles'
		R, V, PPX, PPY, PPZ = load_parts(h,gal_num,r_limit,HaloID,HVX,HVY,HVZ,Z,R_crit200,HVD)

		print '...caustic!'
		x_range,INF_NFWMASS,DIA_NFWMASS,INF_CAUMASS,DIA_CAUMASS,INF_MPROF,INF_NFW,INF_CAU,DIA_MPROF,DIA_NFW,DIA_CAU = P.kernel_caustic_masscalc(R,V,M_crit200,R_crit200,SRAD,ESRAD,HVD,halo_num,bin_range,gal_num,H0,q,r_limit,run_num,use_mems)

		return x_range,INF_NFWMASS,DIA_NFWMASS,INF_CAUMASS,DIA_CAUMASS,INF_MPROF,INF_NFW,INF_CAU,DIA_MPROF,DIA_NFW,DIA_CAU,R,V,HaloID,R_crit200,M_crit200






	
