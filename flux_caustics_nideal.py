'''The caustics module!!!!!!'''
import matplotlib
from matplotlib.pyplot import *
import numpy as np
#import pyfits
import math
#from cosmocalc import cosmocalc    #not for use on FLUX
import scipy.ndimage as ndi
from scipy.integrate import simps
from scipy.optimize import curve_fit
import astStats
import numpy.ma as ma

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
	    '''
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
            '''
            return (rvalues,vvalues,magvalues)
	
	def Limit_richness(self,r,v,mag,r200,N=100,mems=None,limit=None):
            "Limit the sample to only the N brightest galaxies"
            rcheck = r[np.argsort(mag.T[2])]
            vcheck = v[np.argsort(mag.T[2])]
            magcheck = mag[np.argsort(mag.T[2])]
            if limit is None:
                if mems is not None:
                    memcheck = mems[np.argsort(mag.T[2])]
                    return(rcheck[:N],vcheck[:N],magcheck[:N],memcheck[:N])
                else:	# Taking M brightest galaxies such that there are N galaxies within r200. M>=N
                    end = N-1
                    self.check = 0
                    while self.check < N:
                        end += 1
                        self.check = (rcheck[:end])[np.where(rcheck[:end]<r200)].size
                        if end > 250: break
                    return (rcheck[:end],vcheck[:end],magcheck[:end])
            if limit is not None:
                set = np.arange(0,rcheck[magcheck.T[2]<limit].size,1)
                np.random.shuffle(set)
                if rcheck[magcheck.T[2]<limit].size == 0:
                    return (np.array([]),np.array([]),np.array([]))
                else:
                    return (rcheck[magcheck.T[2]<limit][set][:N],vcheck[magcheck.T[2]<limit][set][:N],magcheck[magcheck.T[2]<limit][set][:N])
                
	
	def Selection(self,r,v,mag,r200):
            "Selects galaxies based on a percent criterion in radial bins"
            mag_cut = np.max(mag)-0.01
            N = 100
            while N > 100:
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
                perc = [1.0,0.9,0.8,0.8,0.8]
                #perc = [.6,.6,.6,.8,1.0]
                #perc = [.5,.5,.5,.5,.5]
                #perc = [1.0,1.0,1.0,1.0,1.0]
			
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
	
	def masscalc(self,ri,A,r200,halom,vdisp,density=None,density_tot=None,beta=None,conc=None,dstyle=None):
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
                self.g_b = np.zeros(self.r2.size) + (3-2*0.2)/(1-0.2)
            else:
                self.g_b = ((3.0-2.0*beta)/(1.0-beta))[ri>=0]
            self.f_beta = 0.5*((self.r2/r200)*self.conc)**2/((1+((self.r2/r200)*self.conc))**2*np.log(1+((self.r2/r200)*self.conc)))*self.g_b
            #self.f_beta[0] = 0
            if dstyle==None:
                for i in range(self.A2.size-1):
                    i += 1    
                    self.sum[i] = np.trapz(self.f_beta[1:i+1]*(self.A2[1:i+1]*1000)**2,(self.r2[1:i+1])*self.kmMpc*1000)
                    #self.sum[i] = np.trapz((A2[:i+1]*1000)**2,(r2[:i+1])*kmMpc*1000)
            if dstyle==True:
                for i in range(self.A2.size-1):
                    i += 1
                    self.sum[i] = np.trapz(0.5*(self.A2[1:i+1]*1000)**2,(self.r2[1:i+1])*self.kmMpc*1000)
                
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
            return self.massprof,self.f_beta*(self.A2/3.1e19)**2/(4.5e-48)
	
	def plotcluster(self,r,v,rm,vm,x,y,dens,A,bin,rich_lim,maxv,dens_grad,dens_inf,beta,name):
                g_b = (3.0-2.0*beta)/(1.0-beta)
                g_bc = (3.0-2.0*beta[2])/(1.0-beta[2])
		s = figure()
		ax = s.add_subplot(111)
		ax.pcolormesh(x,y,dens.T)
		#ax.plot(x,np.abs(A)/np.sqrt(g_b),lw=3,c='red') #caustic from main program (probably nfw mphi)
                #ax.plot(x,np.abs(self.Ar_final),c='red',lw=1,ls='--') #original mphi caustic
                ax.plot(x,np.abs(self.Ar_finalD),c='orange',lw=1,ls='--') #diaferio caustic
		#ax.plot(rm,vm,'ro',ls='None',alpha=0.5,markersize=12)
		ax.plot(r,v,'o',color='orange',ls='None',alpha=.8,markersize=5)
                #ax.plot(realbins,np.sqrt(-2*real_potential)*3.08e19/np.sqrt(g_bc),c='red',lw=2)
                #ax.plot(x,np.sqrt(-2*pot_NFW)*3.08e19/np.sqrt(g_bc),c='green',lw=2)
                #ax.plot(x[2:],-np.sqrt(-2*potential[2:])*3.08e19,c='orange',lw=2)
		ax.set_xlim(0,2.5)
		ax.set_ylim(0,3500)
                xlabel('r (Mpc)',fontsize='large')
                ylabel('relative velocity to halo (km/s)',fontsize='large')
		#ax.set_xlabel('r (Mpc)',fontsize='large')
		#s.savefig('/nfs/christoq_ls/giffordw/flux_figs/caustics/nideal/'+str(bin-1)+name+'_gal.png')
		#s.savefig('figures/'+str(bin-1)+'.'+str(rich_lim)+'n.png')
		#show()
		close()

                s = figure()
		ax = s.add_subplot(111)
		ax.pcolormesh(x,y,dens_grad.T)
		ax.plot(x,np.abs(A)/np.sqrt(g_b),lw=3,c='red')
		ax.plot(x,-np.abs(A),lw=3,c='red')
                ax.plot(x,np.abs(self.Ar_final),c='red',lw=1,ls='--')
                ax.plot(x,np.abs(self.Ar_finalD),c='green',lw=1,ls='--')
		#ax.plot(rm,vm,'ro',ls='None',alpha=0.5,markersize=12)
		ax.plot(r,v,'go',ls='None',alpha=.8,markersize=5)
                #ax.plot(x[2:],np.sqrt(-2*potential[2:])*3.08e19/np.sqrt(g_b[2:]),c='orange',lw=2)
		ax.set_xlim(0,2.5)
		ax.set_ylim(0,3500)
		#ax.set_xlabel('r (Mpc)',fontsize='large')
		#s.savefig('/nfs/christoq_ls/giffordw/flux_figs/caustics/nideal/'+str(bin-1)+name+'_grad.png')
		#s.savefig('figures/'+str(bin-1)+'.'+str(rich_lim)+'n.png')
		#show()
		close()

                s = figure()
		ax = s.add_subplot(111)
		ax.pcolormesh(x,y,dens_inf.T)
		ax.plot(x,np.abs(A)/np.sqrt(g_b),lw=3,c='red')
		ax.plot(x,-np.abs(A),lw=3,c='red')
                ax.plot(x,np.abs(self.Ar_final),c='red',lw=1,ls='--')
                ax.plot(x,np.abs(self.Ar_finalD),c='green',lw=1,ls='--')
		#ax.plot(rm,vm,'ro',ls='None',alpha=0.5,markersize=12)
		ax.plot(r,v,'go',ls='None',alpha=.8,markersize=5)
                #ax.plot(x[2:],np.sqrt(-2*potential[2:])*3.08e19/np.sqrt(g_b[2:]),c='orange',lw=2)
		ax.set_xlim(0,2.5)
		ax.set_ylim(0,3500)
		#ax.set_xlabel('r (Mpc)',fontsize='large')
		#s.savefig('/nfs/christoq_ls/giffordw/flux_figs/caustics/nideal/'+str(bin-1)+name+'_inf.png')
		#s.savefig('figures/'+str(bin-1)+'.'+str(rich_lim)+'n.png')
		#show()
		close()
	
	def gaussian_kernel(self,xvalues,yvalues,r200,normalization=100,scale=10,res=200,adj=20,see=False):
	        yres = 220
		#x_scale = (xvalues-np.min(xvalues))/np.max(xvalues-np.min(xvalues))*res
		#y_scale = ((yvalues-np.min(yvalues))/(normalization*scale))/np.max(xvalues-np.min(xvalues))*res
		self.x_scale = xvalues/6.0*res
		self.y_scale = ((yvalues+5000)/(normalization*scale))/(10000.0/(normalization*scale))*yres

		#img = np.zeros((int(np.max(x_scale))+1,int(np.max(y_scale))+1))
		img = np.zeros((res+1,yres+1))
		#x_range = np.linspace(np.min(xvalues),np.max(xvalues),int(np.max(x_scale))+1)
		#y_range = np.linspace(np.min(yvalues),np.max(yvalues),int(np.max(y_scale))+1)
		x_range = np.linspace(0,6,res+1)
		y_range = np.linspace(-5000,5000,yres+1) 
	
		for j in range(xvalues.size):
			img[self.x_scale[j],self.y_scale[j]] += 1
		#pcolormesh(img.T)
		#find ksize
		#xval = xvalues[np.where((xvalues<3) & (yvalues<2000) & (yvalues > -2000))]
		#yval = yvalues[np.where((xvalues<3) & (yvalues<2000) & (yvalues > -2000))]
		#x_scale2 = (xval-np.min(xval))/np.max(xval-np.min(xval))*res
		#y_scale2 = ((yval-np.min(yval))/(normalization*scale))/np.max(xval-np.min(xval))*res
		#xksize = 3.12/(xvalues.size)**(1.0/6.0)*((np.var(x_scale))/2.0)**0.5/adj
		#yksize = 3.12/(xvalues.size)**(1.0/6.0)*((np.var(y_scale))/2.0)**0.5/adj

		self.ksize = 3.12/(xvalues.size)**(1/6.0)*((np.var(self.x_scale[xvalues<r200])+np.var(self.y_scale[xvalues<r200]))/2.0)**0.5/adj
                self.ksize_x = (4.0/(3.0*xvalues.size))**(1/5.0)*np.std(self.x_scale[xvalues<r200])
                self.ksize_y = (4.0/(3.0*yvalues.size))**(1/5.0)*np.std(self.y_scale[xvalues<r200])
                if self.ksize < 3.5:
                    self.ksize = 3.5
		#ksize = 6.77588630223
		#print 'kernel size',ksize
		#img = ndi.uniform_filter(img, (self.ksize,self.ksize))#,mode='reflect')
                img = ndi.gaussian_filter(img, (self.ksize_y,self.ksize_x))#,mode='reflect')
                img_grad = ndi.gaussian_gradient_magnitude(img, (self.ksize_y,self.ksize_x))
                img_inf = ndi.gaussian_gradient_magnitude(ndi.gaussian_gradient_magnitude(img, (self.ksize_y,self.ksize_x)), (self.ksize_y,self.ksize_x))

#		if see == True:
#			s = figure()
#			ax = s.add_subplot(111)
#			ax.pcolormesh(x_range,y_range,img.T)
#			show()
		return (x_range,y_range,img,np.abs(img_grad),np.abs(img_inf))

	def level_search(self,r,v,rmems,vmems,mags,ri,vi,Zi,norm,r200,rlimit,maxv,beta,halo_srad,halo_esrad,use_vdisp=None,use_mems=None,bin=None):
            kappaguess = np.max(Zi)   #first thing is to guess at the level
            levels = np.linspace(0.00001,kappaguess,100)[::-1] #create levels (kappas) to try out
            fitting_radii = np.where((ri>=r200/3.0) & (ri<=r200))
            g_b = (3-2*beta)/(1-beta)
            
            #Here are the conditions if using membership, a fed value, or estimated membership to get vdispersion
            if use_mems is not None: 
                vvar = self.membervdisp(rmems,vmems,vi,ri,r200)
                print 'Members velocity dispersion: %.2f'%(np.sqrt(vvar))
            if use_vdisp is not None:
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
            self.vesc = np.zeros(levels.size)
            Ar_final_opt = np.zeros((levels.size,ri[np.where((ri<r200) & (ri>=0))].size))
            print 'begin loop'
            for i in range(self.vesc.size): # find the escape velocity for all level (kappa) guesses
                self.vesc[i],Ar_final_opt[i] = self.findvesc(levels[i],ri,vi,Zi,norm,r200,g_b)
                #if i > 0 and np.abs(vesc[i]-4*vvar) > np.abs(vesc[i-1]-4*vvar):
                #	break
            skr = (self.vesc-4.0*vvar)**2
            #print np.vstack((vesc,np.zeros(vesc.size)+4*3*vvar,np.zeros(vesc.size)+vesc_true)).T
            try:
                level_elem = np.where(skr == np.min(skr[np.isfinite(skr)]))[0][0]
                #print 'done with loop'
                self.Ar_finalD = np.zeros(ri.size)
                level_final = levels[level_elem]
                #plot(levels,skr)
                #ylim(1e13,5e14)
                #axvline(level_final,c='green')
                #savefig('minimize.png')
                #close()
                print level_final
                for k in range(self.Ar_finalD.size):
                    self.Ar_finalD[k] = self.findAofr(level_final,Zi[k],vi)
                    if k != 0:
                        #Ar_final[k] = norm*restrict_gradient(np.abs(Ar_final[k-1])/norm,np.abs(Ar_final[k])/norm,ri[k-1],ri[k])
                        self.Ar_finalD[k] = self.restrict_gradient2(np.abs(self.Ar_finalD[k-1]),np.abs(self.Ar_finalD[k]),ri[k-1],ri[k])
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
                self.Ar_finalD = np.zeros(ri.size)

            min_func = lambda x,d0: np.sqrt(2*4*np.pi*4.5e-48*d0*(halo_srad)**2*np.log(1+x/halo_srad)/(x/halo_srad))*3.08e19
            v0 = np.array([1e15])
            self.out = curve_fit(min_func,ri[fitting_radii],self.Ar_finalD[fitting_radii],v0[:],maxfev=2000)
            dens_fit = self.out[0][0]
            try:
                self.e_dens = np.sqrt(self.out[1][0][0])
            except:
                self.e_dens = 1e15
            srad_fit = halo_srad
            self.e_srad = halo_esrad
            print 'Used Table Scale Radius'
            vesc_fit = np.sqrt(2*4*np.pi*4.5e-48*dens_fit*(srad_fit)**2*np.log(1+ri/srad_fit)/(ri/srad_fit))*3.08e19
            '''
            for j in range(levels.size):
                plot(ri[np.where((ri<r200) & (ri>=0))],np.abs(Ar_final_opt[j]))
            plot(r,np.abs(v),'k.')
            plot(ri,np.abs(self.Ar_finalD),lw=2,c='green')
            plot(ri,np.abs(self.Ar_final_trueD),lw=2,c='red')
            plot(ri[2:],np.sqrt(-2*potential[2:])*3.08e19,c='blue',lw=2)
            axhline(np.sqrt(4*3*vvar),c='black',ls='--',lw=2)
            axhline(np.sqrt(vesc[level_elem]),c='green',ls='--',lw=2)
            axhline(np.sqrt(vesc_true),c='red',ls='--',lw=2)
            ylim(0,4500)
            savefig('/nfs/christoq_ls/giffordw/flux_figs/surfacetests/'+str(bin-1)+'.png')
            close()
            '''
            print 'done with caustic search'
            return (vesc_fit,3.5*(vvar)**.5,dens_fit,self.e_dens,srad_fit,self.e_srad)

        def level_search2(self,r,v,rmems,vmems,mags,ri,vi,Zi,Zi_inf,norm,r200,rlimit,maxv,beta,halo_srad,halo_esrad,use_vdisp=None,use_mems=None,bin=None):
            kappaguess = np.max(Zi)   #first thing is to guess at the level
            fitting_radii = np.where((ri>=r200/3.0) & (ri<=r200))
            c_guess = np.array([halo_srad])#np.linspace(1.0,12.0,100)
            density_guess = np.linspace(1e13,5e16,1000)
            self.levels = np.linspace(0.00001,kappaguess,150)[::-1] #create levels (kappas) to try out
            g_b = (3-2*beta)/(1-beta)
            #Here are the conditions if using membership, a fed value, or estimated membership to get vdispersion
            if use_mems is not None: 
                vvar = self.membervdisp(rmems,vmems,vi,ri,r200)
                print 'Members velocity dispersion: %.2f'%(np.sqrt(vvar))
            if use_vdisp is not None:
                vvar = use_vdisp**2
                print 'The velocity dispersion was fed from main program: %.2f'%(use_vdisp)
            else:
                vvar = self.findvdisp(r,v,vi,ri,r200,maxv)  #find the variance of the galaxy data
                print 'my clipped vdisp: %.2f'%(np.sqrt(vvar))
            self.Ar_final_opt = np.zeros((self.levels.size,ri[np.where((ri<r200) & (ri>=0))].size))
            self.inf_vals = np.zeros((self.levels.size,ri[np.where((ri<r200) & (ri>=0))].size))
#            s = figure()
#            ax = s.add_subplot(111)
            for i in range(self.levels.size): # find the escape velocity for all level (kappa) guesses
                self.Ar_final_opt[i],self.inf_vals[i] = self.findvesc2(self.levels[i],ri,vi,Zi,Zi_inf,norm,r200)
#                ax.plot(ri[np.where((ri<r200) & (ri>=0))],np.abs(self.Ar_final_opt[i]),c='black',alpha=0.4)
            self.inf_avg = np.average(self.inf_vals.T[fitting_radii],axis=0) #average inflection along each caustic surface
            self.Ar_avg = np.average((self.Ar_final_opt.T[ri<r200]).T,axis=1)
            tryfit = np.polyfit(self.levels,self.inf_avg,7)
            self.infyvals = tryfit[0]*self.levels**7+tryfit[1]*self.levels**6+tryfit[2]*self.levels**5+tryfit[3]*self.levels**4+tryfit[4]*self.levels**3+tryfit[5]*self.levels**2+tryfit[6]*self.levels+tryfit[7]
            '''
            s2 = figure()
            ax2 = s2.add_subplot(111)
            ax2.plot(self.Ar_avg,self.infyvals)
            ax2.plot(self.Ar_avg,self.inf_avg,'k.')
            savefig('/nfs/christoq_ls/giffordw/hi.png')
            close()
            '''
            self.inf_std = np.std(self.inf_vals.T[fitting_radii],axis=0) #std of inflection along each caustic surface
            #self.level_elem = (self.levels[Ar_avg>np.sqrt(vvar)])[np.where(self.inf_avg[Ar_avg>np.sqrt(vvar)] == np.max(self.inf_avg[Ar_avg>np.sqrt(vvar)]))]
            self.level_elem = self.levels[np.where(self.inf_avg == np.max(self.inf_avg))][0]
            #low_zone = np.where((np.average(np.abs(self.Ar_final_opt),axis=1)>np.max(v)/2.0) & (np.average(np.abs(self.Ar_final_opt),axis=1)<np.max(v)))
            high_zone = np.where((np.average(np.abs(self.Ar_final_opt),axis=1)>np.max(v)/2.0))
            #level_elem_low = self.levels[low_zone][np.where(self.inf_avg[low_zone] == np.min(self.inf_avg[low_zone]))][-1]
            #level_elem_high = self.levels[high_zone][np.where(self.inf_avg[high_zone] == np.max(self.inf_avg[high_zone]))][-1]
            try:
                self.level_elem_high = (self.levels[1:-1][np.where((self.infyvals[1:-1]>self.infyvals[2:])&(self.infyvals[1:-1]>self.infyvals[:-2]))])[-1]
            except IndexError:
                self.level_elem_high = self.levels[0]
            self.Ar_final_high = np.zeros(ri.size)
            #self.Ar_final_low = np.zeros(ri.size)
            for i in range(ri.size):
                self.Ar_final_high[i] = self.findAofr(self.level_elem_high,Zi[i],vi)
                #self.Ar_final_low[i] = self.findAofr(level_elem_low,Zi[i],vi)
                if i > 0:
                    self.Ar_final_high[i] = self.restrict_gradient2(np.abs(self.Ar_final_high[i-1]),np.abs(self.Ar_final_high[i]),ri[i-1],ri[i])
                    #self.Ar_final_low[i] = self.restrict_gradient2(np.abs(self.Ar_final_low[i-1]),np.abs(self.Ar_final_low[i]),ri[i-1],ri[i])
            #Ar_final = self.Ar_final_opt[np.where(self.inf_avg == np.max(self.inf_avg))][0]
            #self.Ar_final = (self.Ar_final_high+self.Ar_final_low)/2.0
            self.Ar_final = self.Ar_final_high
            '''
            #########################################
            # If you want a fixed halo concentration#
            #########################################
            min_func = lambda x,d0: np.sqrt(2*4*np.pi*4.5e-48*d0*(halo_srad)**2*np.log(1+x/halo_srad)/(x/halo_srad))*3.08e19
            v0 = np.array([1e15])
            self.out = curve_fit(min_func,ri[fitting_radii],self.Ar_final[fitting_radii],v0[:],maxfev=2000)
            dens_fit = self.out[0][0]
            try:
                self.e_dens = np.sqrt(self.out[1][0][0])
            except:
                self.e_dens = 1e15
            srad_fit = halo_srad
            self.e_srad = halo_esrad
            print 'Used Table Scale Radius'
            '''
            #########################################################
            #If you want to vary both normalization and scale radius#
            #########################################################
            
            print 'BEGIN CHISQ NFW FIT'
            #try because it may not find the correct parameters
            try:
                min_func = lambda x,d0,s0: np.sqrt(2*4*np.pi*4.5e-48*d0*(s0)**2*np.log(1+x/s0)/(x/s0))*3.08e19
                v0 = np.array([1e15,0.5])
                self.out = curve_fit(min_func,ri[fitting_radii],self.Ar_final[fitting_radii]*np.sqrt(g_b[fitting_radii]),v0[:],maxfev=2000)
                dens_fit = self.out[0][0]
                self.e_dens = np.sqrt(self.out[1][0][0])
                srad_fit = self.out[0][1]
                self.e_srad = np.sqrt(self.out[1][1][1])
                print 'Fit Scale Radius'
            
            except RuntimeError:
                srad_fit = r200/2.0
                self.e_srad = srad_fit/2.0
                min_func = lambda x,d0: np.sqrt(2*4*np.pi*4.5e-48*d0*(srad_fit)**2*np.log(1+x/srad_fit)/(x/srad_fit))*3.08e19
                v0 = np.array([1e15])
                self.out = curve_fit(min_func,ri[fitting_radii],self.Ar_final[fitting_radii]*np.sqrt(g_b[fitting_radii]),v0[:],maxfev=2000)
                dens_fit = self.out[0][0]
                try:
                    self.e_dens = np.sqrt(self.out[1][0][0])
                except:
                    self.e_dens = 1e15
            
            if srad_fit > 1.0:
                srad_fit = r200/2.0
                self.e_srad = srad_fit/2.0
                min_func = lambda x,d0: np.sqrt(2*4*np.pi*4.5e-48*d0*(srad_fit)**2*np.log(1+x/srad_fit)/(x/srad_fit))*3.08e19
                v0 = np.array([1e15])
                self.out = curve_fit(min_func,ri[fitting_radii],self.Ar_final[fitting_radii]*np.sqrt(g_b[fitting_radii]),v0[:],maxfev=2000)
                dens_fit = self.out[0][0]
                try:
                    self.e_dens = np.sqrt(self.out[1][0][0])
                except:
                    self.e_dens = 1e15
            if srad_fit < 0.05:
                srad_fit = 0.05
                self.e_srad = srad_fit/2.0
                min_func = lambda x,d0: np.sqrt(2*4*np.pi*4.5e-48*d0*(srad_fit)**2*np.log(1+x/srad_fit)/(x/srad_fit))*3.08e19
                v0 = np.array([1e15])
                self.out = curve_fit(min_func,ri[fitting_radii],self.Ar_final[fitting_radii]*np.sqrt(g_b[fitting_radii]),v0[:],maxfev=2000)
                dens_fit = self.out[0][0]
                try:
                    self.e_dens = np.sqrt(self.out[1][0][0])
                except:
                    self.e_dens = 1e15 
            ri2 = ri[1:]
            #print 'integrate profile', 4*np.pi*dens_fit*(r200/conc_fit)**3*(np.log(1+conc_fit)-conc_fit/(1+conc_fit))
            #print 'in level search',dens_fit,conc_fit
            vesc_fit = np.sqrt(2*4*np.pi*4.5e-48*dens_fit*(srad_fit)**2*np.log(1+ri/srad_fit)/(ri/srad_fit))*3.08e19
#            ax.plot(ri,np.abs(self.Ar_final),c='red',lw=2)
#            ax.plot(ri,vesc_fit,c='green',lw=2)
#            ax.plot(r,v,'k.')
            #pcolormesh(ri,vi,Zi_inf.T)
#            ax.set_ylim(0,3500)
            #savefig('/nfs/christoq_ls/giffordw/flux_figs/surfacetests/nideal/'+str(bin-1)+'.png')
#            close()
            return vesc_fit,3.5*(vvar)**.5,dens_fit,self.e_dens,srad_fit,self.e_srad

	def findvesc(self,level,ri,vi,Zi,norm,r200,g_b):
		'''Calculate vesc^2 by first calculating the integrals in Diaf 99 which is not labeled but in 
		between Eqn 18 and 19'''
		useri = ri[np.where((ri<r200) & (ri>=0))] #look only inside r200
		Ar = np.zeros(useri.size)
		phir = np.zeros(useri.size)
		#loop through each dr and find the caustic amplitude for the given level (kappa) passed to this function
		for i in range(useri.size):
			Ar[i] = self.findAofr(level,Zi[np.where((ri<r200) & (ri>=0))][i],vi)
			if i > -1:  #to fix the fact that the first row of Zi is 'nan'
				#The Serra paper also restricts the gradient when the ln gradient is > 2. We use > 3
				#Ar[i] = norm*restrict_gradient(np.abs(Ar[i-1])/norm,np.abs(Ar[i])/norm,useri[i-1],useri[i])
				Ar[i] = self.restrict_gradient2(np.abs(Ar[i-1]),np.abs(Ar[i]),useri[i-1],useri[i])
				philimit = np.abs(Ar[i]) #phi integral limits
				phir[i] = self.findphir(Zi[i][np.where((vi<philimit) & (vi>-philimit))],vi[np.where((vi<philimit) & (vi>-philimit))])
		return (np.trapz(Ar**2*phir,useri)/np.trapz(phir,useri),Ar)
                #return (np.trapz(g_b[ri<r200]*Ar**2*phir,useri)/np.trapz(phir,useri),Ar)

        def findvesc2(self,level,ri,vi,Zi,Zi_inf,norm,r200):
            useri = ri[np.where((ri<r200) & (ri>=0))] #look only inside r200
            Ar = np.zeros(useri.size)
            inf_val = np.zeros(useri.size)
            for i in range(useri.size):
                Ar[i] = self.findAofr(level,Zi[np.where((ri<r200) & (ri>=0))][i],vi)
                if i >0:
                    Ar[i] = self.restrict_gradient2(np.abs(Ar[i-1]),np.abs(Ar[i]),useri[i-1],useri[i])
                inf_val[i] = Zi_inf[i][np.where(np.abs(vi-Ar[i]) == np.min(np.abs(vi-Ar[i])))][0]
            return Ar,inf_val

	def findAofr(self,level,Zi,vgridvals):
		"Finds the velocity where kappa is"
		dens0 = Zi[np.where(vgridvals>=0)][0]
		if dens0 >= level:
			maxdens = 0.0 #v value we are centering on
			highvalues = Zi[np.where(vgridvals >= maxdens)] #density values above the center v value maxdens
			lowvalues = Zi[np.where(vgridvals < maxdens)] #density values below the center v value maxdens
			highv = vgridvals[np.where(vgridvals >= maxdens)] #v values above the center v value maxdens
			lowv = vgridvals[np.where(vgridvals < maxdens)] #v values below the center v value maxdens
			highslot = self.identifyslot(highvalues,level) #identify the velocity
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
                            if i != 0:
				slot = i-1
				break
                            else:
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

        def twogapper(self,data,gap=300.0):
            avg_r = np.average(data[:,0])
            data = data[np.argsort(data[:,1])]
            #split into 4 quadrants
            databin1 = data[np.where((data[:,0]<avg_r)&(data[:,1]>0))]
            self.databin1 = data[np.where((data[:,0]<avg_r)&(data[:,1]>0))]
            databin2 = data[np.where((data[:,0]<avg_r)&(data[:,1]<=0))]
            databin3 = data[np.where((data[:,0]>=avg_r)&(data[:,1]>0))]
            databin4 = data[np.where((data[:,0]>=avg_r)&(data[:,1]<=0))]
            #calculate each bin's max v gap
            gap1 = databin1[:,1][1:]-databin1[:,1][:-1]
            gap2 = databin2[:,1][1:]-databin2[:,1][:-1]
            gap3 = databin3[:,1][1:]-databin3[:,1][:-1]
            gap4 = databin4[:,1][1:]-databin4[:,1][:-1]

            try:
                if np.max(gap1) >= gap: vgap1 = np.where(gap1 >= gap)[0][0]
                else: vgap1 = 9999999
            except ValueError:
                vgap1 = None
            try:
                if np.max(gap2) >= gap: vgap2 = np.where(gap2 >= gap)[0][-1]
                else: vgap2 = 0
            except ValueError:
                vgap2 = None
            try:
                if np.max(gap3) >= gap: vgap3 = np.where(gap3 >= gap)[0][0]
                else: vgap3 = 9999999
            except ValueError:
                vgap3 = None
            try:
                if np.max(gap4) >= gap: vgap4 = np.where(gap4 >= gap)[0][-1]
                else: vgap4 = 0
            except ValueError:
                vgap4 = None
            print vgap1,vgap2,vgap3,vgap4
            #remove galaxies
            if vgap1 is not None: 
                final_data = databin1[:vgap1+1]
            if vgap2 is not None:
                try: 
                    final_data = np.append(final_data,databin2[vgap2:],axis=0)
                except:
                    final_data = databin2[vgap2:]
            if vgap3 is not None:
                try:
                    final_data = np.append(final_data,databin3[:vgap3+1],axis=0)
                except:
                    final_data = databin3[:vgap3+1]
            if vgap4 is not None:
                try:
                    final_data = np.append(final_data,databin4[vgap4:],axis=0)
                except:
                    final_data = databin4[vgap4:]
            return(final_data)

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

            thetavec = np.arccos(z/np.sqrt(x**2.0+y**2.0+z**2.0))
            phivec = np.arctan(y/x)
            vrad = vx*np.sin(thetavec)*np.cos(phivec)+vy*np.sin(thetavec)*np.sin(phivec)+vz*np.cos(thetavec)
            vtheta = vx*np.cos(thetavec)*np.cos(phivec)+vy*np.cos(thetavec)*np.sin(phivec)-vz*np.sin(thetavec)
            vphi = -vx*np.sin(phivec)+vy*np.cos(phivec)
            rvec = np.sqrt(x**2.0+y**2.0+z**2.0)
            self.beta = np.zeros(radii.size)
            self.beta -= 999.0
            for i in range(radii.size-1):
                i += 1
                w = np.where((rvec>radii[i-1]) & (rvec<=radii[i]))
                if w[0].size >= 20:
                    self.beta[i] = 1.0 - (astStats.biweightScale(vtheta[w],9.0)**2.0 + astStats.biweightScale(vphi[w],9.0)**2.0)/(2.0*astStats.biweightScale(vrad[w],9.0)**2.0)
            #fit = np.polyfit(radii[np.where((self.beta>-5))],self.beta[np.where((self.beta>-5))],6)
            #self.yfit = fit[0]*radii**6.0 + fit[1]*radii**5.0 + fit[2]*radii**4.0 + fit[3]*radii**3.0 + fit[4]*radii**2.0 + fit[5]*radii + fit[6]
            return self.beta
