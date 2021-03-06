'''The caustics module!!!!!!'''
from matplotlib.pyplot import *
import numpy as np
#import pyfits
import math
#from cosmocalc import cosmocalc    #not for use on FLUX
import scipy.ndimage as ndi
import astStats

class caustic:
	def angulardistance(self,clus_z,H=100.0):
		'''Finds the angular diameter distance for an array of cluster center redshifts. Cannot use on FLUX.
		Instead, use angular distance file precalculated and upload.'''
		#print 'begin anglular distance calc'
		try:
			ang_d = [cosmocalc(z,H0=H)['DA_Mpc'] for z in clus_z]  #in Mpc
		except TypeError:
			ang_d = cosmocalc(clus_z,H0=H)['DA_Mpc']
		#for i in range(len(clus_z)): ang_d = c*np.array(clus_z)/H
		return ang_d

	def findangle(self,ra,dec,clus_RA,clus_DEC):
		'''This function takes 4 arrays and an index. The first input is the index of the center you are 
		looking at. The second/third input are ra/dec arrays for a list of objects you want to know how 
		far away they are from the center. The fourth/fifth input are the ra/dec arrays for the centers.
		angle is returned in radians'''
		zsep = np.sin(clus_DEC*np.pi/180.0)*np.sin(np.array(dec)*np.pi/180.0)
		xysep = np.cos(clus_DEC*np.pi/180.0)*np.cos(np.array(dec)*math.pi/180.0)*np.cos(np.pi/180.0*(clus_RA-np.array(ra)))
		angle = np.arccos(zsep+xysep)
		return angle
	
	def set_sample(self,r,v,mags,rlimit,vlimit=3500,H0=72.0,gal_mem=None):
            ''' The gal_mem argument allows for members to be returned so you have an extra returning value'''
            rvalues = np.array(r)[np.where((r < rlimit) & (v < vlimit) & (v > -vlimit))]
            vvalues = np.array(v)[np.where((r < rlimit) & (v < vlimit) & (v > -vlimit))]
            magvalues = np.array(mags)[np.where((r < rlimit) & (v < vlimit) & (v > -vlimit))]
		
            try: #This fixes the 0 velocity location based on the N galaxies inside our limits
                vfix = astStats.biweightLocation(vvalues[np.where((rvalues<0.5) & (vvalues>-vlimit) & (vvalues<vlimit))],6.0)
                vvalues = vvalues - vfix
            except: #Exception is caught because astStats needs a certain number of galaxies to work
                vfix = np.average(vvalues[np.where((rvalues<0.5) & (vvalues>-3500) & (vvalues<3500))])
                vvalues = vvalues - vfix
            rvalues1 = np.array(rvalues)[np.where((rvalues < rlimit) & (vvalues < vlimit) & (vvalues > -vlimit))]
            vvalues1 = np.array(vvalues)[np.where((rvalues < rlimit) & (vvalues < vlimit) & (vvalues > -vlimit))]
            magvalues1 = np.array(magvalues)[np.where((rvalues < rlimit) & (vvalues < vlimit) & (vvalues > -vlimit))]

            if gal_mem is not None: #If gal_mem is called, it will return an extra value
                memvalues = np.array(gal_mem)[np.where((r < rlimit) & (v < vlimit) & (v > -vlimit))]
                memvalues1 = np.array(memvalues)[np.where((rvalues < rlimit) & (vvalues < vlimit) & (vvalues > -vlimit))]
                return (rvalues1,vvalues1,magvalues1,memvalues1)
            else:
                return (rvalues1,vvalues1,magvalues1)
	
	def Limit_richness(self,r,v,mag,r200,N=100,mems=None):
            "Limit the sample to only the N brightest galaxies"
            rcheck = r[np.argsort(mag.T[2])]
            vcheck = v[np.argsort(mag.T[2])]
            magcheck = mag[np.argsort(mag.T[2])]
            if mems is not None:
                memcheck = mems[np.argsort(mag.T[2])]
                return(rcheck[:N],vcheck[:N],magcheck[:N],memcheck[:N])
            else:
                return (rcheck[:N],vcheck[:N],magcheck[:N])
	
	def Selection(self,r,v,mag,r200):
            "Selects galaxies based on a percent criterion in radial bins"
            mag_cut = np.max(mag)-0.01
            N = 100
            while N > 50:
                vin = v[np.where(mag<mag_cut)]
                rin = r[np.where(mag<mag_cut)]
                magin = mag[np.where(mag<mag_cut)]
                (n,bins) = np.histogram(rin,bins=5)
                d = np.digitize(rin,bins[:-1])
		    
                ########################
                #Choose percent profile#
                ########################
                #example cored percent profile [.6,.6,.6,.8,1.0]
                #example winged percent profile [1.0,0.8,0.6,0.6,0.6]
                #example uniform percent profile [.7,.7,.7,.7,.7]
                #perc = [1.0,0.8,0.6,0.6,0.6]
                #perc = [.6,.6,.6,.8,1.0]
                #perc = [.5,.5,.5,.5,.5]
                perc = [1.0,1.0,1.0,1.0,1.0]
			
                r = np.array([])
                v = np.array([])
                mag = np.array([])
                for i in range(n.size):
                    i+=1 #because d is in terms of bin number not element
                    r = np.append(rin[np.where(d==i)][np.argsort(magin[np.where(d==i)])[:np.ceil(rin[np.where(d==i)].size*perc[i-1])]],r)
                    v = np.append(vin[np.where(d==i)][np.argsort(magin[np.where(d==i)])[:np.ceil(rin[np.where(d==i)].size*perc[i-1])]],v)
                    mag = np.append(magin[np.where(d==i)][np.argsort(magin[np.where(d==i)])[:np.ceil(rin[np.where(d==i)].size*perc[i-1])]],mag)
            N = v.size
            mag_cut -= 0.01
            print N
            return (r,v)
	
	def MissCenter(self,HaloRa,HaloDec,z):
		"Offset the cluster center by a given amount"
		self.offset = 1.5
		HaloRa += (self.offset/np.sqrt(2))/angulardistance(z)*180.0/np.pi
		HaloDec += (self.offset/np.sqrt(2))/angulardistance(z)*180.0/np.pi
		return HaloRa,HaloDec
		
	def Verror(self,v,sigv):
		"Apply an extra error to the velocity"
		rand = np.random.normal(0,sigv,v.size)
		return v + rand
	
	def masscalc(self,ri,A,r200,halom,vdisp,density=None,density_tot=None,beta=None,conc=None):
            "Calculate the mass profile"
            self.G = 6.67E-11
            self.per = 8.6e8
            self.solmass = 1.98892e30
            self.r2 = ri[ri>=0]
            self.A2 = A[ri>=0]
            self.kmMpc = 3.08568025e19
            self.sum = np.zeros(self.A2.size)
            #if beta == None:
            #conc = 5*(halom/1e14)**-0.1
            if conc == None:
                self.conc = 4.0*(vdisp/700.0)**(-0.306)
            else:
                self.conc = conc
            print 'concentration = ', self.conc
            if beta == None:
                self.g_b = np.zeros(self.r2.size) + (3-2*0.3)/(1-0.3)
            else:
                self.g_b = (3.0-2.0*beta)/(1.0-beta)
            self.f_beta = 0.5*((self.r2/r200)*self.conc)**2/((1+((self.r2/r200)*self.conc))**2*np.log(1+((self.r2/r200)*self.conc)))*self.g_b
            self.f_beta[0] = 0
            for i in range(self.A2.size-1):
                i += 1    
                self.sum[i] = np.trapz(self.f_beta[:i+1]*(self.A2[:i+1]*1000)**2,(self.r2[:i+1])*self.kmMpc*1000)
                #self.sum[i] = np.trapz((A2[:i+1]*1000)**2,(r2[:i+1])*kmMpc*1000)
            self.massprof = self.sum/(self.G*self.solmass)
            '''
            else:
                self.g_b = (3.0-2.0*beta)/(1.0-beta)
                self.pot = np.zeros(self.r2.size)
                for i in range(self.r2.size-1):
                    i += 1
                    #self.pot[i] = self.G*(density_tot[i]*4.0/3.0*np.pi*((self.r2[i]*self.kmMpc*1000)**2.0))
                    self.pot[i] = self.G/(self.r2[i]*self.kmMpc*1000)*np.trapz(4*np.pi*(density_tot[1:i+1]*self.solmass/(self.kmMpc*1000)**3)*(self.r2[1:i+1]*self.kmMpc*1000)**2,self.r2[1:i+1]*self.kmMpc*1000)+4*np.pi*np.trapz(density_tot[i:]*self.solmass/(self.kmMpc*1000)**3*self.r2[i:]*self.kmMpc*1000,self.r2[i:]*self.kmMpc*1000)
                    self.sum[i] = np.trapz(self.G*2.0*np.pi*(self.A2[1:i+1]*1000)**2.0*self.g_b[1:i+1]*(density_tot[1:i+1]*self.solmass/(self.kmMpc*1000)**3.0)*(self.r2[1:i+1]*self.kmMpc*1000)**2.0/self.pot[1:i+1],self.r2[1:i+1]*self.kmMpc*1000)
                    #self.sum[i] = np.trapz(3*self.g_b[1:i+1]*(self.A2[1:i+1]*1000)**2/2.0,self.r2[1:i+1]*self.kmMpc*1000)
            self.massprof = self.sum/(self.solmass*self.G)
            '''
            return self.massprof
	
	def plotcluster(self,r,v,rm,vm,x,y,dens,A,bin,rich_lim,maxv):
		s = figure()
		ax = s.add_subplot(111)
		ax.pcolormesh(x,y,dens.T)
		ax.plot(x,np.abs(A),lw=3,c='red')
		ax.plot(x,-np.abs(A),lw=3,c='red')
		#ax.plot(rm,vm,'ro',ls='None',alpha=0.5,markersize=12)
		ax.plot(r,v,'go',ls='None',alpha=.8,markersize=6)
		ax.set_xlim(0,5)
		ax.set_ylim(-3500,3500)
		#ax.set_xlabel('r (Mpc)',fontsize='large')
		#ax.set_ylabel(r'v$_{pec}$ (km/s)',fontsize='large')
		s.savefig('figures/caustics/'+str(rich_lim)+'n/'+str(bin)+'.'+str(rich_lim)+'n.png')
		#show()
		close()
	
	def gaussian_kernel(self,xvalues,yvalues,r200,normalization=100,scale=10,res=200,adj=20,see=False):
		yres = 200
		#x_scale = (xvalues-np.min(xvalues))/np.max(xvalues-np.min(xvalues))*res
		#y_scale = ((yvalues-np.min(yvalues))/(normalization*scale))/np.max(xvalues-np.min(xvalues))*res
		x_scale = xvalues/6.0*res
		y_scale = ((yvalues+4000)/(normalization*scale))/(8000.0/(normalization*scale))*yres

		#img = np.zeros((int(np.max(x_scale))+1,int(np.max(y_scale))+1))
		img = np.zeros((res+1,yres+1))
		#x_range = np.linspace(np.min(xvalues),np.max(xvalues),int(np.max(x_scale))+1)
		#y_range = np.linspace(np.min(yvalues),np.max(yvalues),int(np.max(y_scale))+1)
		x_range = np.linspace(0,6,res+1)
		y_range = np.linspace(-4000,4000,yres+1) 
	
		for j in range(xvalues.size):
			img[x_scale[j],y_scale[j]] += 1
		#pcolormesh(img.T)
		#find ksize
		#xval = xvalues[np.where((xvalues<3) & (yvalues<2000) & (yvalues > -2000))]
		#yval = yvalues[np.where((xvalues<3) & (yvalues<2000) & (yvalues > -2000))]
		#x_scale2 = (xval-np.min(xval))/np.max(xval-np.min(xval))*res
		#y_scale2 = ((yval-np.min(yval))/(normalization*scale))/np.max(xval-np.min(xval))*res
		#xksize = 3.12/(xvalues.size)**(1.0/6.0)*((np.var(x_scale))/2.0)**0.5/adj
		#yksize = 3.12/(xvalues.size)**(1.0/6.0)*((np.var(y_scale))/2.0)**0.5/adj
	
		ksize = 3.12/(xvalues.size)**(1/6.0)*((np.var(x_scale[xvalues<r200])+np.var(y_scale[xvalues<r200]))/2.0)**0.5/adj
		#ksize = 6.77588630223
		#print 'kernel size',ksize
		img = ndi.gaussian_filter(img, (ksize,ksize))#,mode='reflect')
		
		if see == True:
			s = figure()
			ax = s.add_subplot(111)
			ax.pcolormesh(x_range,y_range,img.T)
			show()
		return (x_range,y_range,img)

	def level_search(self,r,v,rmems,vmems,mags,ri,vi,Zi,norm,r200,rlimit,maxv,use_vdisp=False,use_mems=False):
            kappaguess = np.max(Zi)   #first thing is to guess at the level
            levels = np.linspace(0.00001,kappaguess,100)[::-1] #create levels (kappas) to try out
            
            #Here are the conditions if using membership, a fed value, or estimated membership to get vdispersion
            if use_mems: 
                vvar = self.membervdisp(rmems,vmems,vi,ri,r200)
                print 'Members velocity dispersion: %.2f'%(np.sqrt(vvar))
            elif use_vdisp:
                vvar = use_vdisp**2
                print 'The velocity dispersion was fed from main program: %.2f'%(use_vdisp)
            else:
                vvar = self.findvdisp(r,v,vi,ri,r200,maxv)  #find the variance of the galaxy data
                print 'my clipped vdisp: %.2f'%(np.sqrt(vvar))
                '''
                try:
                    vvar = (astStats.biweightClipped(v,9.0,3.0)['biweightScale'])**2.0
                    print 'Estimated membership velocity dispersion: %.2f'%(np.sqrt(vvar))
                except:
                    vvar = np.var(v)
                    print 'Estimated membership velocity dispersion: %.2f'%(np.sqrt(vvar))
                '''
            #vvar = findvdisp2(r,v,mags,vi,ri,r200,maxv)  #use radial bins to do sigma clipping and find variance (DON'T USE FOR NOW)
            #vvar = findvdisp3(r,v,mags,r200,maxv) #use red sequence to identify members (DON'T USE FOR NOW)
            vesc = np.zeros(levels.size)
            print 'begin loop'
            for i in range(vesc.size): # find the escape velocity for all level (kappa) guesses
                vesc[i] = self.findvesc(levels[i],ri,vi,Zi,norm,r200)
                #if i > 0 and np.abs(vesc[i]-4*vvar) > np.abs(vesc[i-1]-4*vvar):
                #	break
            skr = (vesc-4*vvar)**2
            try:
                level_elem = np.where(skr == np.min(skr[np.isfinite(skr)]))[0][0]
                #print 'done with loop'
                Ar_final = np.zeros(ri.size)
                level_final = levels[level_elem]
                print level_final
                for k in range(Ar_final.size):
                    Ar_final[k] = self.findAofr(level_final,Zi[k],vi)
                    if k != 0:
                        #Ar_final[k] = norm*restrict_gradient(np.abs(Ar_final[k-1])/norm,np.abs(Ar_final[k])/norm,ri[k-1],ri[k])
                        Ar_final[k] = self.restrict_gradient2(np.abs(Ar_final[k-1]),np.abs(Ar_final[k]),ri[k-1],ri[k])
                    #only if you want
                    '''
                    A_kappa = np.zeros((vesc.size,Ar_final.size))
                    for i in range(levels.size):
                        for k in range(Ar_final.size):
                            A_kappa[i][k] = findAofr(levels[i],Zi[k],vi)
                            if k != 0:
                                A_kappa[i][k] = norm*restrict_gradient(np.abs(A_kappa[i][k-1])/norm,np.abs(A_kappa[i][k])/norm,ri[k-1],ri[k])
                    area = np.trapz(np.abs(A_kappa),ri)
                    print 'final level', level_final
                    s = plt.figure()
                    ax = s.add_subplot(111)
                    ax.plot(levels,area)
                    ax.axvline(x=level_final)
                    #plt.xlabel(r"$\kappa$")
                    #plt.ylabel("area under caustic")
                    plt.show()
                    #plt.close()
                    #s.savefig('/n/Pictor1/giffordw/research/pysim/virialcond.eps')
                    '''
            except ValueError:  #This exception occurs if skr is entirely NAN. A flag should be raised for this in the output table
                Ar_final = np.zeros(ri.size)
            print 'done with caustic search'
            return Ar_final,3.5*(vvar)**.5

	def findvesc(self,level,ri,vi,Zi,norm,r200):
		'''Calculate vesc^2 by first calculating the integrals in Diaf 99 which is not labeled but in 
		between Eqn 18 and 19'''
		useri = ri[np.where((ri<r200) & (ri>=0))] #look only inside r200
		Ar = np.zeros(useri.size)
		phir = np.zeros(useri.size)
		#loop through each dr and find the caustic amplitude for the given level (kappa) passed to this function
		for i in range(useri.size):
			Ar[i] = self.findAofr(level,Zi[np.where((ri<r200) & (ri>=0))][i],vi)
			if i != 0:  #to fix the fact that the first row of Zi is 'nan'
				#The Serra paper also restricts the gradient when the ln gradient is > 2. We use > 3
				#Ar[i] = norm*restrict_gradient(np.abs(Ar[i-1])/norm,np.abs(Ar[i])/norm,useri[i-1],useri[i])
				Ar[i] = self.restrict_gradient2(np.abs(Ar[i-1]),np.abs(Ar[i]),useri[i-1],useri[i])
				philimit = np.abs(Ar[i]) #phi integral limits
				phir[i] = self.findphir(Zi[i][np.where((vi<philimit) & (vi>-philimit))],vi[np.where((vi<philimit) & (vi>-philimit))])
		return np.trapz(Ar[1:]**2*phir[1:],useri[1:])/np.trapz(phir[1:],useri[1:])

	def findAofr(self,level,Zi,vgridvals):
		"Finds the velocity where kappa is"
		dens0 = Zi[np.where(vgridvals>=0)][0]
		if dens0 >= level:
			maxdens = 0.0
			highvalues = Zi[np.where(vgridvals >= maxdens)]
			lowvalues = Zi[np.where(vgridvals < maxdens)]
			highv = vgridvals[np.where(vgridvals >= maxdens)]
			lowv = vgridvals[np.where(vgridvals < maxdens)]
			highslot = self.identifyslot(highvalues,level)
			flip_lowslot = self.identifyslot(lowvalues[::-1],level)
			lowslot = lowvalues.size - flip_lowslot
			if len(lowv) == 0 or len(highv) == 0: #probably all zeros
				highamp = lowamp = 0
				return highamp
			if highslot == highv.size:
				highamp = highv[-1]
			if lowslot ==0:
				lowamp = lowv[0]
			if highslot == 0 or lowslot == lowv.size:
				highamp = lowamp = 0
			if highslot != 0 and highslot != highv.size:
				highamp = highv[highslot]-(highv[highslot]-highv[highslot-1])*(1-(highvalues[highslot-1]-level)/(highvalues[highslot-1]-highvalues[highslot]))
			if lowslot != 0 and lowslot != lowv.size:
				lowamp = lowv[lowslot-1]-(lowv[lowslot-1]-lowv[lowslot])*(1-(lowvalues[lowslot]-level)/(lowvalues[lowslot]-lowvalues[lowslot-1]))
			if np.abs(highamp) >= np.abs(lowamp):
				return lowamp
			if np.abs(highamp) < np.abs(lowamp):
				return highamp
		else: return 0 #no maximum density exists

	def restrict_gradient(self,pastA,newA,pastr,newr):
		if pastA <= newA:
			if (newA-pastA)/(newr-pastr) > 3.0:
				#print newA,pastA
				dr = newr-pastr
				return pastA + 0.25*dr
			else: return newA
		if pastA > newA:
			if (newA-pastA)/(newr-pastr) < -3.0:
				#print newA,pastA
				dr = newr-pastr
				return pastA - 0.25*dr
			else: return newA

	def restrict_gradient2(self,pastA,newA,pastr,newr):
		"It is necessary to restrict the gradient the caustic can change at in order to be physical"
		if pastA <= newA:
			if (np.log(newA)-np.log(pastA))/(np.log(newr)-np.log(pastr)) > 3.0:
				#print newA,pastA
				dr = np.log(newr)-np.log(pastr)
				return np.exp(np.log(pastA) + 2*dr)
			else: return newA
		if pastA > newA:
			if (np.log(newA)-np.log(pastA))/(np.log(newr)-np.log(pastr)) < -3.0 and pastA != 0:
				#print newA,pastA
				dr = np.log(newr)-np.log(pastr)
				return np.exp(np.log(pastA) - 2*dr)
			else: return newA

	def findphir(self,shortZi,shortvi):
		short2Zi = np.ma.masked_array(shortZi)
		#print 'test',shortvi[np.ma.where(np.ma.getmaskarray(short2Zi)==False)]
		vi = shortvi[np.ma.where(np.ma.getmaskarray(short2Zi)==False)]
		Zi = short2Zi[np.ma.where(np.ma.getmaskarray(short2Zi)==False)]
		
		vi = vi[np.isfinite(Zi)]
		Zi = Zi[np.isfinite(Zi)]
		x = np.trapz(Zi.compressed(),vi)
		return x

	def identifyslot(self,dvals,level):
		'''This function takes the density values for a given r grid value either above or below
		the v grid value that corresponds to the maximum density at the r slice and returns the indici
		where the level finally falls below the given level. Density values should be in order
		starting with the corresponding value to the v value closest to the maximum and working toward
		the edges (high to low density in general).'''
		slot = dvals.size - 1
		for i in range(dvals.size):
			if level >= dvals[i]:
				slot = i
				break
		return slot

	def findvdisp(self,r,v,vi,ri,r200,maxv):
            #print 'average r', np.average(r)
            avgr = r200
            #dispvals = v[np.where((r>np.average(r)-.4) & (r<np.average(r)+.4) & (v<2000) & (v>-2000))]
            for i in range(6):
                v2 = v[np.where((r<avgr) & (v<maxv) & (v>-maxv))]
                r2 = r[np.where((r<avgr) & (v<maxv) & (v>-maxv))]
                stv = 3.5 * np.std(v2)
                print '3.5 sigma of v = ', stv
                v = v2[np.where((v2 > -stv) & (v2 < stv))]
                r = r2[np.where((v2 > -stv) & (v2 < stv))]
            if v.size > 15.0:
                vstd = astStats.biweightScale(v,9.0)
                vvar = (astStats.biweightScale(v,9.0))**2
            else:
                vstd = np.std(v)
                vvar = np.var(v)
            #print 'standard dev of zone= ',vstd
            return (np.sqrt(vvar))**2

	def findvdisp2(self,r,v,mags,vi,ri,r200,maxv):
		"do radial sigma clipping to find vdisp"
		binedge = np.arange(0,r200+1,0.3)
		rin = r[np.where((r<r200) & (v<maxv) & (v>-maxv))]
		vin = v[np.where((r<r200) & (v<maxv) & (v>-maxv))]
		vfinal = np.array([])
		for i in range(binedge.size-1):
			i += 1
			x = rin[np.where((rin>binedge[i-1]) & (rin<binedge[i]))]
			y = vin[np.where((rin>binedge[i-1]) & (rin<binedge[i]))]
			for k in range(6):
				y2 = y
				x2 = x
				stv = 3.5 * np.std(y2)
				#print '3.5 sigma of v = ',stv
				y = y2[np.where((y2 > -stv) & (y2 < stv))]
				x = x2[np.where((y2 > -stv) & (y2 < stv))]
			vstd2 = np.std(y)
			vvar2 = np.var(y)
			print 'standard dev of zone %i = %f' % (i,vstd2)
			vfinal = np.append(y[np.where((y<vvar2) & (y>-vvar2))],vfinal)
		return np.var(vfinal)

	def findvdisp3(self,r,v,mags,r200,maxv):
		"use red sequence to find members"
		binedge = np.arange(0,r200+1,0.3)
		rin = r
		vin = v
		colin = mags.T[1] - mags.T[2]
		avg_c = np.average(colin)
		vfinal = np.array([])
		for i in range(binedge.size-1):
			i += 1
			x = rin[np.where((rin>binedge[i-1]) & (rin<binedge[i]))]
			y = vin[np.where((rin>binedge[i-1]) & (rin<binedge[i]))]
			c = colin[np.where((rin>binedge[i-1]) & (rin<binedge[i]))]
			for k in range(6):
				y2 = y
				x2 = x
				c2 = c
				stv = 3.5 * np.std(y2)
				y = y2[np.where((y2 > -stv) & (y2 < stv) | ((c2<avg_c+0.04) & (c2>avg_c-0.04)))]
				x = x2[np.where((y2 > -stv) & (y2 < stv) | ((c2<avg_c+0.04) & (c2>avg_c-0.04)))]
				c = c2[np.where((y2 > -stv) & (y2 < stv) | ((c2<avg_c+0.04) & (c2>avg_c-0.04)))]
			vstd2 = np.std(y)
			vvar2 = np.var(y)
			print 'standard dev of zone %i = %f' % (i,vstd2)
			vfinal = np.append(y[np.where((y<vvar2) & (y>-vvar2))],vfinal)
		return np.var(vfinal)

	def findvdisp4(self,r,v,r200,maxv):
		"shifting gapper method"
		k = False
		b = 6
		while k == False:
			b -= 1
			(n,bins) = np.histogram(r,bins=b)
			k = np.all([n>15])
		print 'bin sizes', n
		d = np.digitize(r,bins[:-1])
		v_final = np.array([])
		r_final = np.array([])
		for i in range(n.size):
			velocities_p = np.sort(v[np.where((d==i+1) & (v>0))])
			radius_p = (r[np.where((d==i+1) & (v>0))])[np.argsort(v[np.where((d==i+1) & (v>0))])]
			velocities_n = np.sort(v[np.where((d==i+1) & (v<0))])[::-1]
			radius_n = (r[np.where((d==i+1) & (v<0))])[np.argsort(v[np.where((d==i+1) & (v<0))])[::-1]]
			dv_p = velocities_p[1:] - velocities_p[:-1]
			dv_n = velocities_n[:-1] - velocities_n[1:]
			for j in range(dv_p.size):
				if dv_p[j] >= 1000.0:
					v_final = np.append(v_final,velocities_p[:j+1])
					r_final = np.append(r_final,radius_p[:j+1])
					break
			for j in range(dv_n.size):
				if dv_n[j] >= 1000.0:
					v_final = np.append(v_final,velocities_n[:j+1])
					r_final = np.append(r_final,radius_n[:j+1])
					break
		try:
			vvar = (astStats.biweightScale(v,9.0))**2
		except:
			vvar = np.var(v)
		return vvar
	
	
	def membervdisp(self,r,v,vi,ri,r200):
		"This function is for the ideal scenario that you know which galaxies are members"
		#print 'standard dev of zone= ',np.std(v[np.where((r<r200))])# & (v>-2000) & (v < 2000))])
		#return np.var(v)
		#return np.var(v[np.where((r<r200) & (v>-2000) & (v < 2000))])
		try:
			vvar = (astStats.biweightScale(v,9.0))**2.0
		except:
			vvar = np.var(v)
		return vvar

        def densityprofile(self,x,y,z,halox,haloy,haloz,radii):
            density = np.zeros(radii.size)
            density_tot = np.zeros(radii.size)
            r = np.sqrt((x-halox)**2+(y-haloy)**2+(z-haloz)**2)
            r_cut = r[np.where(r<=np.max(radii))]
            for i in range(radii.size-1):
                i+=1
                density[i] = 8.6e8*r_cut[np.where((r_cut<radii[i]) & (r_cut>radii[i-1]))].size/(4.0/3.0*np.pi*radii[i]**3.0-4.0/3.0*np.pi*radii[i-1]**3.0)
                density_tot[i] = 8.6e8*r_cut[np.where(r_cut<radii[i])].size/(4.0/3.0*np.pi*radii[i]**3.0)
            return (density,density_tot)

        def betaprofile(self,x,y,z,vx,vy,vz,halox,haloy,haloz,halovx,halovy,halovz,radii,rlimit):
            #go to cluster reference frame
            x = x-halox
            y = y-haloy
            z = z-haloz
            #correct for cluster proper motion
            vx = vx-halovx
            vy = vy-halovy
            vz = vz-halovz

            thetavec = np.arccos(vz/np.sqrt(vx**2.0+vy**2.0+vz**2.0))
            phivec = np.arctan(vx/vy)
            vrad = vx*np.sin(thetavec)*np.cos(phivec)+vy*np.sin(thetavec)*np.sin(phivec)+vz*np.cos(thetavec)
            vtheta = vx*np.cos(thetavec)*np.cos(phivec)+vy*np.cos(thetavec)*np.sin(phivec)+vz*np.sin(thetavec)
            vphi = -vx*np.sin(phivec)+vy*np.cos(phivec)
            rvec = np.sqrt(x**2.0+y**2.0+z**2.0)
            self.beta = np.zeros(radii.size)
            self.beta -= 999.0
            for i in range(radii.size-1):
                i += 1
                w = np.where((rvec>radii[i-1]) & (rvec<=radii[i]))
                if w[0].size >= 20:
                    self.beta[i] = 1.0 - (astStats.biweightScale(vtheta[w],9.0)**2.0 + astStats.biweightScale(vphi[w],9.0)**2.0)/(2.0*astStats.biweightScale(vrad[w],9.0)**2.0)
            fit = np.polyfit(radii[np.where((self.beta>-5))],self.beta[np.where((self.beta>-5))],7)
            self.yfit = fit[0]*radii**7.0 + fit[1]*radii**6.0 + fit[2]*radii**5.0 + fit[3]*radii**4.0 + fit[4]*radii**3.0 + fit[5]*radii**2.0 + fit[6]*radii + fit[7]
            return self.yfit
