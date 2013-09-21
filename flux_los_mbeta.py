'''finding masses of Millenium clusters'''
import numpy as np
from flux_caustics_ideal import *
import pyfits
import sys
import astStats
from matplotlib.pyplot import *
from scipy import weave
from scipy.weave import converters
import cosmolopy.distance as cd

class get_data:
    def get_gals(self,ID,H0):
        "Specify the directory/filename(s) of the galaxies in the simulation"
        galaxy_xpos,galaxy_ypos,galaxy_zpos,galaxy_vx,galaxy_vy,galaxy_vz,galaxy_umag,galaxy_gmag,galaxy_rmag,galaxy_imag,galaxy_zmag = np.loadtxt('minimill/minidata/halo_'+ID+'_d35_all_galaxies.csv',dtype='float',delimiter=',',usecols=(1,2,3,4,5,6,7,8,9,10,11),unpack=True)
        galaxy_id = np.loadtxt('minimill/minidata/halo_'+ID+'_d35_all_galaxies.csv',dtype='string',delimiter=',',usecols=(0,),unpack=True)
        return (galaxy_id,galaxy_umag,galaxy_gmag,galaxy_rmag,galaxy_imag,galaxy_zmag,galaxy_xpos/(H0/100.0),galaxy_ypos/(H0/100.0),galaxy_zpos/(H0/100.0),galaxy_vx,galaxy_vy,galaxy_vz)

    def get_galsbig(self,ID,H0,Z):
        galaxy_z,galaxy_umag,galaxy_gmag,galaxy_rmag,galaxy_imag,galaxy_zmag,galaxy_xpos,galaxy_ypos,galaxy_zpos,galaxy_vx,galaxy_vy,galaxy_vz = np.loadtxt('data/lowz_data2_2/'+ID+'.galaxies.tab',dtype='float',usecols=(3,4,5,6,7,8,9,10,11,12,13,14),unpack=True)
        galaxy_id = np.loadtxt('data/lowz_data2_2/'+ID+'.galaxies.tab',dtype='string',usecols=(0,),unpack=True)
        galaxy_xpos,galaxy_ypos,galaxy_zpos = galaxy_xpos/(1+Z),galaxy_ypos/(1+Z),galaxy_zpos/(1+Z)
        return (galaxy_id,galaxy_z,galaxy_umag,galaxy_gmag,galaxy_rmag,galaxy_imag,galaxy_zmag,galaxy_xpos/(H0/100.0),galaxy_ypos/(H0/100.0),galaxy_zpos/(H0/100.0),galaxy_vx,galaxy_vy,galaxy_vz)

    def get_galsbig2(self,ID,H0):
        fileID = np.loadtxt('/nfs/christoq_ls/MILLENNIUM/particles/cmiller.csv',dtype='string',delimiter=',',skiprows=1,usecols=(0,),unpack=True)
        Nid = np.where(ID==fileID)[0][0]
        fits = pyfits.open('/nfs/christoq_ls/MILLENNIUM/particles/t_'+str(Nid)+'_cmiller_guo.fits')
        data = fits[1].data
        return (data.field('HALOID'),data.field('uDUST'),data.field('gDUST'),data.field('rDUST'),data.field('iDUST'),data.field('zDUST'),data.field('x')/(H0/100.0),data.field('y')/(H0/100.0),data.field('z')/(H0/100.0),data.field('velX'),data.field('velY'),data.field('velZ'),data.field('vvir'))
        

    def get_clusters(self,H0):
        "Specify the directory/filename of the halos in the simulation"
        M200,clus_vdisp,r200,cluster_xpos,cluster_ypos,cluster_zpos,cluster_vx,cluster_vy,cluster_vz = np.loadtxt('minimill/halos.csv',dtype='float',delimiter=',',usecols=(4,16,7,12,13,14,17,18,19),unpack=True)
        haloid,subhaloid = np.loadtxt('minimill/halos.csv',dtype='string',delimiter=',',usecols=(1,2),unpack=True)
        #ang_distances = np.loadtxt('data/ang_distances_lowz2.tab',dtype='float',usecols=(0,),unpack=True)
        return (haloid,M200,subhaloid,cluster_xpos/(H0/100.0),cluster_ypos/(H0/100.0),cluster_zpos/(H0/100.0),cluster_vx,cluster_vy,cluster_vz,r200,clus_vdisp)

    def get_clustersbig(self,H0):
        "Specify the directory/filename of the halos in the simulation"
        r200,M200,clus_vdisp,clus_xpos,clus_ypos,clus_zpos,clus_vx,clus_vy,clus_vz = np.loadtxt('data/biglosclusters.csv',dtype='float',delimiter=',',usecols=(5,6,8,9,10,11,12,13,14),unpack=True)
        haloid,subhaloid = np.loadtxt('data/biglosclusters.csv',dtype='string',delimiter=',',usecols=(0,1),unpack=True)
        scale_rad,e_scale_rad,scale_dens,e_scale_dens,r200,r500,M200crit,M500crit,clus_vdisp,clus_vdisp500,snaps,red,hofz = np.loadtxt('/nfs/christoq_ls/giffordw/Millenium/files/nfwprofs/Millbig_concentrations.phys_phys.csv',dtype='float',delimiter=',',usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13),unpack=True)
        return (haloid,M200,M200crit,subhaloid,clus_xpos/(1+red),clus_ypos/(1+red),clus_zpos/(1+red),clus_vx,clus_vy,clus_vz,r200,r500,clus_vdisp,scale_rad,e_scale_rad,scale_dens,e_scale_dens,clus_vdisp500,snaps,red,hofz)

    def Pick_pos(self,halop):
        "Picks a random position for the observer a given distance away from the center of the chosen halo position"
        x = np.random.uniform(-1,1)
        y = np.random.uniform(-1,1)
        z = np.random.uniform(-1,1)
        unit = np.array([x,y,z])/(x**2+y**2+z**2)**(.5)
        # move the position randomly 50Mpc away
        return halop+30*unit

if __name__ == "__main__":
    #define the classes
    G = get_data()
    C = caustic()

    #set constants
    H0 = 100.0
    q = 10.0
    c = 300000.0
    cosmo = {'omega_M_0':0.3,'omega_lambda_0':0.7,'h':H0/100.0}
    cosmo = cd.set_omega_k_0(cosmo)

    #load halo data from file by calling the get_clusters function
    print 'GETTING CLUSTERS'
    haloid,M200,M200crit,subhaloid,clus_xpos,clus_ypos,clus_zpos,clus_vx,clus_vy,clus_vz,r200,r500,clus_vdisp,clus_srad,clus_esrad,clus_scale,clus_escale,clus_vdisp500,snaps,red,hofz = G.get_clustersbig(H0)
    #the millimill vdisp isn't crit in table so I used particles to calculate
    #vdisp_crit = np.loadtxt('big_vdisp_crit.csv',dtype='float',delimiter=',',usecols=(1,),unpack=True)

    halo_p = np.array([clus_xpos,clus_ypos,clus_zpos])
    halo_v = np.array([clus_vx,clus_vy,clus_vz])
    N = haloid.size
    
    #Both False: Use Dan's sigma clipping -- 
    #Use_vdisp = True: Use the table vdisp -- 
    #use_mems = True: Use known members for vdisp (outdated)
    use_mems = False
    use_vdisp = True
    
    line_num = 1 #This specifies how many lines of sight for each cluster
    rich_lim = 1000000 #This specifies how many galaxies will be used to estimate the caustic surface

    #All my data files are sorted by M200_crit, so here is the sorting block
    HaloID = haloid[np.argsort(M200)[::-1]][:N]
    HaloR200 = r200[np.argsort(M200)[::-1]][:N]/(H0/100.0)
    HaloR500 = r500[np.argsort(M200)[::-1]][:N]/(H0/100.0)
    HaloM200 = M200crit[np.argsort(M200)[::-1]][:N]/(H0/100.0)
    HaloVD = clus_vdisp[np.argsort(M200)[::-1]][:N]
    HaloVD500 = clus_vdisp500[np.argsort(M200)[::-1]][:N]
    Halo_PX = clus_xpos[np.argsort(M200)[::-1]][:N]
    Halo_PY = clus_ypos[np.argsort(M200)[::-1]][:N]
    Halo_PZ = clus_zpos[np.argsort(M200)[::-1]][:N]
    Halo_VX = clus_vx[np.argsort(M200)[::-1]][:N]
    Halo_VY = clus_vy[np.argsort(M200)[::-1]][:N]
    Halo_VZ = clus_vz[np.argsort(M200)[::-1]][:N]
    Halo_srad = clus_srad[np.argsort(M200)[::-1]][:N]
    Halo_esrad = clus_esrad[np.argsort(M200)[::-1]][:N]
    Halo_scale = clus_scale[np.argsort(M200)[::-1]][:N]
    Halo_escale = clus_escale[np.argsort(M200)[::-1]][:N]
    Halo_snap = snaps[np.argsort(M200)[::-1]][:N]
    Halo_Z = red[np.argsort(M200)[::-1]][:N]
    Halo_hofz = hofz[np.argsort(M200)[::-1]][:N]

    try: #if the user feeds a number at the command line, that halo is used. It should usually be 0.
        ar = np.int(sys.argv[1])
    except IndexError: #if the user does not feed a number, the index is chosen to start with 0
        ar = 0
    gal_vdisp3d = np.zeros(100)
    particle_vdisp3d = np.zeros(100)
    caustic_mass = np.zeros(100)
    for k in range(1):#HaloID.size): #loop over the halos
        i = ar+k
        #i = 0 #halo number to use
        mem_flag = 1 #This is a flag to alerting to (1) if members > 0 and (0) if not. Affects number output.
        HaloZ = Halo_Z[i]
        
        #Get current galaxy info
        print 'GETTING GALAXIES'
        
        gal_haloid,gal_z,gal_umag,gal_gmag,gal_rmag,gal_imag,gal_zmag,gal_xpos,gal_ypos,gal_zpos,gal_vx,gal_vy,gal_vz = G.get_galsbig(HaloID[i],H0,HaloZ)
        gal_mags = np.vstack((gal_umag,gal_gmag,gal_rmag,gal_imag,gal_zmag)).T
        lumdist = cd.luminosity_distance(gal_z,**cosmo)
        gal_abs_rmag = gal_rmag - 5*np.log10(lumdist*1e6/10.0)
        '''
        gal_haloid,gal_umag,gal_gmag,gal_rmag,gal_imag,gal_zmag,gal_xpos,gal_ypos,gal_zpos,gal_vx,gal_vy,gal_vz,gal_vdisp = G.get_galsbig2(HaloID[i],H0)
        gal_mags = np.vstack((gal_umag,gal_gmag,gal_rmag,gal_imag,gal_zmag)).T
        gal_abs_rmag = gal_rmag
        '''
        gal_p = np.array([gal_xpos,gal_ypos,gal_zpos])
        gal_v = np.array([gal_vx,gal_vy,gal_vz])
        gal_mem = np.zeros(gal_umag.size)+1

        #organize the current halo position and velocity
        Halo_P = np.array([Halo_PX[i],Halo_PY[i],Halo_PZ[i]]) #current halo position
        Halo_V = np.array([Halo_VX[i],Halo_VY[i],Halo_VZ[i]]) #current halo velocity
        HVD = HaloVD[i]
        HVD500 = HaloVD500[i]


        gal_radius = np.sqrt((Halo_P[0]-gal_p[0])**2+(Halo_P[1]-gal_p[1])**2+(Halo_P[2]-gal_p[2])**2)
        gal_velocity = np.sqrt((Halo_V[0]-gal_v[0])**2+(Halo_V[1]-gal_v[1])**2+(Halo_V[2]-gal_v[2])**2)
        
        #load wherever you have stored the mass profiles output from flux_los_part
        mprof_xtrue,mprof_true,vesc = np.loadtxt('/nfs/christoq_ls/giffordw/Millenium/files/mprofsbig/'+str(HaloID[i])+'data/'+str(HaloID[i])+'.0.1000000mprof.tab',dtype='float',usecols=(0,2,4),unpack=True)
        potential = (vesc/3.08e19)**2/-2.0

        
        #define which galaxy indicies are members
        #mems = np.where(np.array(gal_haloid) == HaloID[i])

        line_mass = np.zeros(line_num) #empty array for different lines of sight masses
        vdisp_lostrue = np.zeros(line_num)
        vdispersion = np.zeros(line_num) #empty array for different lines of sight vdisp
        rnums = np.zeros(line_num)
        #line_error = np.zeros(line_num) #empty array for the errors above

        for j in range(line_num): #loop over different lines of sight
            #define r200 and limits to select our particles for the caustic estimation
            rand_r200 = HaloR200[i]
            vlimit = 5000
            rlimit = 5#rand_r200*1.25
            '''
            ##############################################
            #All line of sight calculations are done here#
            ##############################################
            new_pos = G.Pick_pos(Halo_P) #Decide where the observer sits
            #derive halo information with new los
            halo_dist = ((Halo_P[0]-new_pos[0])**2 + (Halo_P[1]-new_pos[1])**2 + (Halo_P[2]-new_pos[2])**2)**0.5
            halo_pos_vect = np.array([Halo_P[0]-new_pos[0],Halo_P[1]-new_pos[1],Halo_P[2]-new_pos[2]])/halo_dist
            halo_vlos = np.dot(halo_pos_vect, Halo_V)
            #derive galaxy information with new los
            gal_dist = ((gal_p[0]-new_pos[0])**2 + (gal_p[1]-new_pos[1])**2 + (gal_p[2]-new_pos[2])**2)**0.5
            gal_vlos = np.zeros(gal_dist.size)
            n = gal_dist.size
            gal_pos_vect = np.zeros((3,gal_dist.size))
            code = """
            int u,w;
            for (u=0;u<n;++u){
            for(w=0;w<3;++w){
                gal_pos_vect(w,u) = (gal_p(w,u)-new_pos(w))/gal_dist(u);
                }
            gal_vlos(u) = gal_pos_vect(0,u)*gal_v(0,u)+gal_pos_vect(1,u)*gal_v(1,u)+gal_pos_vect(2,u)*gal_v(2,u);
            }
            """
            fast = weave.inline(code,['gal_pos_vect','n','gal_dist','gal_vlos','gal_v','new_pos','gal_p'],type_converters=converters.blitz,compiler='gcc')
            angles = np.arccos(np.dot(halo_pos_vect,gal_pos_vect))
            '''
            r = gal_radius#angles*halo_dist
            v = gal_velocity#gal_vlos-halo_vlos*np.dot(halo_pos_vect,gal_pos_vect)
            gal_vdisp3d[i] = np.sqrt(astStats.biweightScale(gal_v[0][np.where(gal_radius<=HaloR200[i])]-Halo_V[0],9.0)**2+astStats.biweightScale(gal_v[1][np.where(gal_radius<=HaloR200[i])]-Halo_V[1],9.0)**2+astStats.biweightScale(gal_v[2][np.where(gal_radius<=HaloR200[i])]-Halo_V[2],9.0)**2)/np.sqrt(3)
            #print 'MY VELOCITY OF GALAXIES', gal_vdisp3d[i]
            particle_vdisp3d[i] = HVD*np.sqrt(3)
            gal_rmag_new = gal_abs_rmag# + 5*np.log10(gal_dist*1e6/10.0)
            '''
            rand_r200 = findr200(r,v,gal_rmag_new,angles,gal_lumdist,HaloAD[i],H0)*0.615#*1.1
            vlimit = 3500
            rlimit = rand_r200*1.25
            '''
            #import average beta profile and create a fit. Apply fit to your xvalues later in code
            xbeta,abeta = np.loadtxt('data/average_betaprofile.tab',dtype='float',usecols=(0,1),unpack=True)
            fit = np.polyfit((xbeta*rand_r200)[xbeta<4],abeta[xbeta<4],6)
            

            ###################################
            #End of line of sight calculations#
            ###################################
            
            #limit our particles based on set limits in projected radius, los velocity, and richness (aka: what we would observe)
            if use_mems:
                rvalues,vvalues,magvalues,memvalues = C.set_sample(r,v,gal_mags,rlimit,vlimit,H0,gal_mem=gal_mem)
                rvalues,vvalues,magvalues,memvalues = C.Limit_richness(rvalues,vvalues,magvalues,rlimit,N=rich_lim,mems=memvalues)
            else:
                rvalues,vvalues,magvalues = C.set_sample(r,v,gal_mags,rlimit,vlimit,H0)
                rvalues,vvalues,magvalues = C.Limit_richness(rvalues,vvalues,magvalues,rlimit,N=rich_lim)
            rvalues = np.append(rvalues,rvalues)
            vvalues = np.append(vvalues,-vvalues)
            '''
            if mems[0].size > 0:
    	        rmems,vmems,magmems = C.set_sample(r[mems],v[mems],gal_mags[mems],rlimit,vlimit,H0)
                mem_flag = 1
            else:
    	        rmems,vmems,magmems = (rvalues,vvalues,magvalues)
                mem_flag = 0
                print 'NO MEMBERS'
	    '''	
            #Use radial selection template (Not using now)
            #rmems,vmems = C.Selection(rmems,vmems,magmems.T[2],rand_r200)
            rnums[j] = rvalues.size

            ###########################
            #Kernal Density Estimation#
            ###########################
            res_array = np.array([1])
            img_tot = 0
            img_grad_tot = 0
            img_inf_tot = 0
            for u in range(res_array.size):
                x_range,y_range,img,img_grad,img_inf = C.gaussian_kernel(rvalues,vvalues,rand_r200,normalization=H0,scale=q,res=200,adj = res_array[u],see=False)
                img_tot += img/np.max(img)
                img_grad_tot += img_grad/np.max(img_grad)
                img_inf_tot += img_inf/np.max(img_inf)
            ###########################
            #End of Density Estimation#
            ###########################

            #############
            #Define Beta#
            #############
            '''
            try:
                beta = np.loadtxt('files/mprofsbig/'+str(HaloID[i])+'data/'+str(HaloID[i])+'.beta_nosmooth.tab',dtype='float',usecols=(1,),unpack=True)
            except:
                beta = np.zeros(x_range.size)+0.2
            '''
            beta = fit[0]*x_range**6 + fit[1]*x_range**5 + fit[2]*x_range**4 + fit[3]*x_range**3 + fit[4]*x_range**2 + fit[5]*x_range + fit[6]
            #####################
            #End Beta Definition#
            #####################

            ############################
            #Caustic Surface Estimation#
            ############################
            maxv = rand_r200*H0*np.sqrt(200)+500
            Anew,threesig,dens_norm,e_dens_norm,srad,e_srad = C.level_search2(rvalues,vvalues,rvalues,vvalues,rvalues,x_range,y_range,img_tot,img_inf_tot,H0*q,HaloR200[i],rlimit,maxv,beta,potential,Halo_srad[i],Halo_esrad[i],use_vdisp=HVD,use_mems=use_mems,bin=i+1)
            AnewD,threesigD,dens_normD,e_dens_normD,sradD,e_sradD = C.level_search(rvalues,vvalues,rvalues,vvalues,rvalues,x_range,y_range,img_tot,H0*q,HaloR200[i],rlimit,maxv,beta,potential,Halo_srad[i],Halo_esrad[i],use_vdisp=gal_vdisp3d[i],use_mems=use_mems,bin=i+1) #use for limit richness
            vdispersion[j] = threesig/3.5
            ###########################
            #End of Caustic Estimation#
            ###########################

            M200_int = 4*np.pi*dens_norm*(srad)**3*(np.log(1+rand_r200/srad)-rand_r200/srad/(1+rand_r200/srad))
            M200_intD = 4*np.pi*dens_normD*(sradD)**3*(np.log(1+rand_r200/sradD)-rand_r200/sradD/(1+rand_r200/sradD))
            Gus_mass = 1e15*(gal_vdisp3d[i]/1082.9)**(1/0.3361)
            
            #Calculate the mass profile using the caustic equation. The beta profile is optional.
            massprofile,integrand = C.masscalc(x_range,np.abs(C.Ar_finalD),rand_r200,HaloM200[i],vdispersion[j],beta=beta,conc=rand_r200/Halo_srad[i])
            massprofile2,integrand2 = C.masscalc(x_range,np.abs(C.Ar_final),rand_r200,HaloM200[i],vdispersion[j],beta=beta,conc=rand_r200/Halo_srad[i])

            #Identify the mass within a given radius.
            lookwithin = rand_r200  #how many Mpc to print mass within
            M200_intD2 = massprofile[np.where(x_range[x_range >=0] < lookwithin)[0][-1]]
            M200_int2 = massprofile2[np.where(x_range[x_range >=0] < lookwithin)[0][-1]]
            line_mass[j] = M200_int
            print '   M_phi mass = %e' % (line_mass[j]), 'Halo mass = %e' %(HaloM200[i]),'within %.2f'%(rand_r200)
            print '   NFW C mass = %e' % (M200_intD), 'Halo mass = %e' %(HaloM200[i]),'within %.2f'%(rand_r200)
            print '   Caus. mass = %e' % (M200_intD2), 'Halo mass = %e' %(HaloM200[i]),'within %.2f'%(rand_r200)
            print '     Gus mass = %e' % (Gus_mass), 'Halo mass = %e' %(HaloM200[i]),'within %.2f'%(rand_r200)
            print '     looping ', j+1 #print which line of sight loop we are currently on
            
            ####################
            #Print Mass Profile#
            ####################
            fi = open('/nfs/christoq_ls/giffordw/Millenium/files/mprofsbig/'+str(HaloID[i])+'data/'+str(HaloID[i])+'.'+str(j)+'.'+str(rich_lim)+'mprofgals_3d.tab','w')
            for ii in range(x_range.size):
                fi.write(str(x_range[ii])+'\t')
                fi.write(str(massprofile[ii])+'\t')
                fi.write(str(mprof_true[ii])+'\t')
                fi.write(str(C.Ar_final[ii])+'\t')
                fi.write(str(C.Ar_finalD[ii])+'\t')
                fi.write('\n')
            fi.close()

        
        Mphi_hat_NFW = np.average(line_mass)
        Mphi_hat_orig = M200_int2
        caustic_mass_NFW = M200_intD
        caustic_mass_orig = M200_intD2
        chi_sq_fit = np.sum((Anew[1:]-C.Ar_final[1:])**2/C.Ar_final[1:])
        #Plot the typical caustic diagram for this cluster/los
        C.plotcluster(rvalues,vvalues,rvalues,vvalues,x_range,y_range,img_tot,Anew,i+1,rich_lim,maxv,potential,img_grad_tot,img_inf_tot,'mbeta3d')
        print 'DONE WITH CLUSTER ',i+1
        #output the important data to file. Each halo has a separate file and will need to be combined later.
        f = open('/nfs/christoq_ls/giffordw/Millenium/files/mprofsbig/'+str(HaloID[i])+'data/'+str(HaloID[i])+'.'+str(rich_lim)+'.gal_data.tab','w')
        f.write(str(HaloM200[i])+'\t')
        f.write(str(rand_r200)+'\t')
        f.write(str(rvalues.size)+'\t')
        f.write(str(Mphi_hat_NFW)+'\t')
        f.write(str(Mphi_hat_orig)+'\t')
        f.write(str(caustic_mass_orig)+'\t')
        f.write(str(caustic_mass_NFW)+'\t')
        f.write(str(Gus_mass)+'\t')
        f.write(str(HVD)+'\t')
        f.write(str(dens_norm)+'\t')
        f.write(str(e_dens_norm)+'\t')
        f.write(str(srad)+'\t')
        f.write(str(e_srad)+'\t')
        f.write(str(chi_sq_fit)+'\t')
        f.write('\n')
        f.close()
