'''finding masses of Millenium clusters
This program imports galaxy or particle data from the Millennium simulation catalogs for the 100 halos for which we have 3D particle data. It then projects these galaxies/particles along a random line of sight. The program also has the ability to loop over bootstrap resamplings of the N brightest galaxies chosen to do the analysis. Be sure to check the file outputs at the end. The bootstrap output should only be uncommented if the bootstrapping loop is active.
'''
import numpy as np
from flux_caustics_nideal import *
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

    def get_galsbigGuo(self,ID,H0,Z):
        fileID = np.loadtxt('/n/Christoq1/MILLENNIUM/particles/cmiller.csv',dtype='string',delimiter=',',skiprows=1,usecols=(0,),unpack=True)
        Nid = np.where(ID==fileID)[0][0]
        fits = pyfits.open('/n/Christoq1/MILLENNIUM/particles/t_'+str(Nid)+'_cmiller_guo.fits')
        xpos,ypos,zpos = fits[1].data.field('x')/(1+Z),fits[1].data.field('y')/(1+Z),fits[1].data.field('z')/(1+Z)
        return (fits[1].data.field('HALOID'),fits[1].data.field('uDUST'),fits[1].data.field('gDUST'),fits[1].data.field('rDUST'),fits[1].data.field('iDUST'),fits[1].data.field('zDUST'),xpos/(H0/100.0),ypos/(H0/100.0),zpos/(H0/100.0),fits[1].data.field('velX'),fits[1].data.field('velY'),fits[1].data.field('velZ'),fits[1].data.field('vvir'))

    def get_galsbigDelucia(self,ID,H0,Z):
        fileID = np.loadtxt('/nfs/christoq_ls/MILLENNIUM/particles/cmiller.csv',dtype='string',delimiter=',',skiprows=1,usecols=(0,),unpack=True)
        Nid = np.where(ID==fileID)[0][0]
        fits = pyfits.open('/nfs/christoq_ls/MILLENNIUM/particles/t_'+str(Nid)+'_cmiller_delucia2006a.fits')
        xpos,ypos,zpos = fits[1].data.field('x')/(1+Z),fits[1].data.field('y')/(1+Z),fits[1].data.field('z')/(1+Z)
        return (fits[1].data.field('HALOID'),fits[1].data.field('MAG_BDUST'),fits[1].data.field('MAG_VDUST'),fits[1].data.field('MAG_RDUST'),fits[1].data.field('MAG_IDUST'),fits[1].data.field('MAG_KDUST'),xpos/(H0/100.0),ypos/(H0/100.0),zpos/(H0/100.0),fits[1].data.field('velX'),fits[1].data.field('velY'),fits[1].data.field('velZ'),fits[1].data.field('vvir'))

    def get_galsbigBertone(self,ID,H0,Z):
        fileID = np.loadtxt('/nfs/christoq_ls/MILLENNIUM/particles/cmiller.csv',dtype='string',delimiter=',',skiprows=1,usecols=(0,),unpack=True)
        Nid = np.where(ID==fileID)[0][0]
        fits = pyfits.open('/nfs/christoq_ls/MILLENNIUM/particles/t_'+str(Nid)+'_cmiller_bertone2007a.fits')
        xpos,ypos,zpos = fits[1].data.field('x')/(1+Z),fits[1].data.field('y')/(1+Z),fits[1].data.field('z')/(1+Z)
        return (fits[1].data.field('HALOID'),fits[1].data.field('SDSS_uDUST'),fits[1].data.field('SDSS_gDUST'),fits[1].data.field('SDSS_rDUST'),fits[1].data.field('SDSS_iDUST'),fits[1].data.field('SDSS_zDUST'),xpos/(H0/100.0),ypos/(H0/100.0),zpos/(H0/100.0),fits[1].data.field('velX'),fits[1].data.field('velY'),fits[1].data.field('velZ'),fits[1].data.field('vvir'))

    def get_galsbigBower(self,ID,H0,Z):
        fileID = np.loadtxt('/nfs/christoq_ls/MILLENNIUM/particles/cmiller.csv',dtype='string',delimiter=',',skiprows=1,usecols=(0,),unpack=True)
        Nid = np.where(ID==fileID)[0][0]
        fits = pyfits.open('/nfs/christoq_ls/MILLENNIUM/particles/t_'+str(Nid)+'_cmiller_bower2006a.fits')
        xpos,ypos,zpos = fits[1].data.field('x')/(1+Z),fits[1].data.field('y')/(1+Z),fits[1].data.field('z')/(1+Z)
        return (fits[1].data.field('GalaxyID'),fits[1].data.field('u_SDSS'),fits[1].data.field('g_SDSS'),fits[1].data.field('r_SDSS'),fits[1].data.field('i_SDSS'),fits[1].data.field('z_SDSS'),xpos/(H0/100.0),ypos/(H0/100.0),zpos/(H0/100.0),fits[1].data.field('velX'),fits[1].data.field('velY'),fits[1].data.field('velZ'))

    def get_galsbigSubHalos(self,ID,H0,Z):
        fileID = np.loadtxt('/nfs/christoq_ls/MILLENNIUM/particles/cmiller.csv',dtype='string',delimiter=',',skiprows=1,usecols=(0,),unpack=True)
        Nid = np.where(ID==fileID)[0][0]
        fits = pyfits.open('/nfs/christoq_ls/MILLENNIUM/particles/t_'+str(Nid)+'_cmiller_subhalos.fits')
        xpos,ypos,zpos = fits[1].data.field('x')/(1+Z),fits[1].data.field('y')/(1+Z),fits[1].data.field('z')/(1+Z)
        return (fits[1].data.field('HaloID'),np.random.random(xpos.size),np.random.random(xpos.size),-fits[1].data.field('M_CRIT200'),np.random.random(xpos.size),np.random.random(xpos.size),xpos/(H0/100.0),ypos/(H0/100.0),zpos/(H0/100.0),fits[1].data.field('velX'),fits[1].data.field('velY'),fits[1].data.field('velZ'))
    
    def get_galsbigParticles(self,ID,H0,Z):
        fileID = np.loadtxt('/nfs/christoq_ls/MILLENNIUM/particles/cmiller.csv',dtype='string',delimiter=',',skiprows=1,usecols=(0,),unpack=True)
        snapnum = np.loadtxt('/nfs/christoq_ls/MILLENNIUM/particles/cmiller.csv',dtype='float',delimiter=',',skiprows=1,usecols=(1,),unpack=True)
        Nid = np.where(ID==fileID)[0][0]
        fits = pyfits.open('/nfs/christoq_ls/MILLENNIUM/particles/t'+str(Nid)+'_cmiller.dat.fits')
        part_x = fits[1].data.field('PPX')/(1+Z)
        part_y = fits[1].data.field('PPY')/(1+Z)
        part_z = fits[1].data.field('PPZ')/(1+Z)
        part_vx = fits[1].data.field('VVX')/np.sqrt(1+Z)
        part_vy = fits[1].data.field('VVY')/np.sqrt(1+Z)
        part_vz = fits[1].data.field('VVZ')/np.sqrt(1+Z)
        return (fits[1].data.field('ID'),np.random.random(part_x.size),np.random.random(part_x.size),np.random.random(part_x.size),np.random.random(part_x.size),np.random.random(part_x.size),part_x/(H0/100.0),part_y/(H0/100.0),part_z/(H0/100.0),part_vx,part_vy,part_vz)
        
    
    def get_clusters(self,H0):
        "Specify the directory/filename of the halos in the simulation"
        M200,clus_vdisp,r200,cluster_xpos,cluster_ypos,cluster_zpos,cluster_vx,cluster_vy,cluster_vz = np.loadtxt('minimill/halos.csv',dtype='float',delimiter=',',usecols=(4,16,7,12,13,14,17,18,19),unpack=True)
        haloid,subhaloid = np.loadtxt('minimill/halos.csv',dtype='string',delimiter=',',usecols=(1,2),unpack=True)
        #ang_distances = np.loadtxt('data/ang_distances_lowz2.tab',dtype='float',usecols=(0,),unpack=True)
        return (haloid,M200,subhaloid,cluster_xpos/(H0/100.0),cluster_ypos/(H0/100.0),cluster_zpos/(H0/100.0),cluster_vx,cluster_vy,cluster_vz,r200,clus_vdisp)

    def get_clustersbig(self,H0):
        "Specify the directory/filename of the halos in the simulation"
        M200,clus_vdisp,clus_xpos,clus_ypos,clus_zpos,clus_vx,clus_vy,clus_vz = np.loadtxt('biglosclusters.csv',dtype='float',delimiter=',',usecols=(6,8,9,10,11,12,13,14),unpack=True)
        haloid,subhaloid = np.loadtxt('biglosclusters.csv',dtype='string',delimiter=',',usecols=(0,1),unpack=True)
        scale_rad,e_scale_rad,scale_dens,e_scale_dens,r200,r500,M200crit,M500crit,clus_vdisp,clus_vdisp500,snaps,red,hofz = np.loadtxt('Millbig_concentrations.phys_phys.csv',dtype='float',delimiter=',',usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13),unpack=True)
        return (haloid,M200,M200crit,subhaloid,clus_xpos/(H0/100.0*(1+red)),clus_ypos/(H0/100.0*(1+red)),clus_zpos/(H0/100.0*(1+red)),clus_vx,clus_vy,clus_vz,r200,r500,clus_vdisp,scale_rad,e_scale_rad,scale_dens,e_scale_dens,clus_vdisp500,snaps,red,hofz)

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
    cosmo = {'omega_M_0':0.25,'omega_lambda_0':0.75,'h':H0/100.0}
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
    jacknife = None #set to anything to have jacknife performed. None for line-of sight test    

    if jacknife is not None:
        line_num = 3 #This specifies how many lines of sight for each cluster
    else: line_num = 1		# number of lines of sight per halo
    rich_lim = 100 #This specifies how many galaxies will be used to estimate the caustic surface

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
    semianalytic = 'Guo'
    for k in range(1):#HaloID.size): #loop over the halos
        i = ar+k
        #i = 0 #halo number to use
        mem_flag = 1 #This is a flag to alerting to (1) if members > 0 and (0) if not. Affects number output.
        HaloZ = Halo_Z[i]
        
        #Get current galaxy info
        print 'GETTING GALAXIES'
        '''
        gal_haloid,gal_z,gal_umag,gal_gmag,gal_rmag,gal_imag,gal_zmag,gal_xpos,gal_ypos,gal_zpos,gal_vx,gal_vy,gal_vz = G.get_galsbig(HaloID[i],H0,HaloZ)
        gal_mags = np.vstack((gal_umag,gal_gmag,gal_rmag,gal_imag,gal_zmag)).T
        lumdist = cd.luminosity_distance(gal_z,**cosmo)
        gal_abs_rmag = gal_rmag - 5*np.log10(lumdist*1e6/10.0)
        '''
        if semianalytic == 'Guo': gal_haloid,gal_umag,gal_gmag,gal_rmag,gal_imag,gal_zmag,gal_xpos,gal_ypos,gal_zpos,gal_vx,gal_vy,gal_vz,gal_vdisp = G.get_galsbigGuo(HaloID[i],H0,HaloZ)
        if semianalytic == 'Delucia': gal_haloid,gal_umag,gal_gmag,gal_rmag,gal_imag,gal_zmag,gal_xpos,gal_ypos,gal_zpos,gal_vx,gal_vy,gal_vz,gal_vdisp = G.get_galsbigDelucia(HaloID[i],H0,HaloZ)
        if semianalytic == 'Bertone': gal_haloid,gal_umag,gal_gmag,gal_rmag,gal_imag,gal_zmag,gal_xpos,gal_ypos,gal_zpos,gal_vx,gal_vy,gal_vz,gal_vdisp = G.get_galsbigBertone(HaloID[i],H0,HaloZ)
        if semianalytic == 'Bower': gal_haloid,gal_umag,gal_gmag,gal_rmag,gal_imag,gal_zmag,gal_xpos,gal_ypos,gal_zpos,gal_vx,gal_vy,gal_vz = G.get_galsbigBower(HaloID[i],H0,HaloZ)
        if semianalytic == 'Subhalos': gal_haloid,gal_umag,gal_gmag,gal_rmag,gal_imag,gal_zmag,gal_xpos,gal_ypos,gal_zpos,gal_vx,gal_vy,gal_vz = G.get_galsbigSubHalos(HaloID[i],H0,HaloZ)
        #gal_haloid,gal_umag,gal_gmag,gal_rmag,gal_imag,gal_zmag,gal_xpos,gal_ypos,gal_zpos,gal_vx,gal_vy,gal_vz = G.get_galsbigParticles(HaloID[i],H0,HaloZ)
        gal_mags = np.vstack((gal_umag,gal_gmag,gal_rmag,gal_imag,gal_zmag)).T
        gal_abs_rmag = gal_rmag
        
        gal_p = np.array([gal_xpos,gal_ypos,gal_zpos])
        gal_v = np.array([gal_vx,gal_vy,gal_vz])
        gal_mem = np.zeros(gal_umag.size)+1

        #organize the current halo position and velocity
        Halo_P = np.array([Halo_PX[i],Halo_PY[i],Halo_PZ[i]]) #current halo position
        Halo_V = np.array([Halo_VX[i],Halo_VY[i],Halo_VZ[i]]) #current halo velocity
        HVD = HaloVD[i]
        HVD500 = HaloVD500[i]


        gal_radius = np.sqrt((Halo_P[0]-gal_p[0])**2+(Halo_P[1]-gal_p[1])**2+(Halo_P[2]-gal_p[2])**2)
        #gal_radius = np.sqrt(gal_p[0]**2+gal_p[1]**2+gal_p[2]**2)
        gal_velocity = np.sqrt((Halo_V[0]-gal_v[0])**2+(Halo_V[1]-gal_v[1])**2+(Halo_V[2]-gal_v[2])**2)
        
        #load wherever you have stored the mass profiles output from flux_los_part
        mprof_xtrue,mprof_true,vesc = np.loadtxt('/n/Christoq1/giffordw/Millenium/files/mprofsbig/'+str(HaloID[i])+'data/'+str(HaloID[i])+'.0.1000000mprof.tab',dtype='float',usecols=(0,2,4),unpack=True)
        potential = (vesc/3.08e19)**2/-2.0

        
        #define which galaxy indicies are members
        #mems = np.where(np.array(gal_haloid) == HaloID[i])

        line_mass = np.zeros(line_num) #empty array for different lines of sight masses
        line_vdisp = np.zeros(line_num)
        vdisp_lostrue = np.zeros(line_num)
        vdispersion = np.zeros(line_num) #empty array for different lines of sight vdisp
        rnums = np.zeros(line_num)
        surf_err = np.zeros(line_num)
        #line_error = np.zeros(line_num) #empty array for the errors above
	''' I (Nick) commented this out b/c I don't want to over-write Dan's files
        li = open('/nfs/christoq_ls/giffordw/Millenium/files/mprofsbig/'+str(HaloID[i])+'data/'+str(HaloID[i])+'.'+str(rich_lim)+'.nideal_lines.'+semianalytic+'.tab','w')
        li.write('#Mphi_mass caustic_mass caustic_mass_fbeta gal_vdisp scale_dens e_scale_dens scale_rad e_scale_rad\n')
        '''
        for j in range(line_num): #loop over different lines of sight
            #define r200 and limits to select our particles for the caustic estimation
            rand_r200 = HaloR200[i]
            vlimit = 3500
            rlimit = rand_r200*1.25
            
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
            r = angles*halo_dist
            v = gal_vlos-halo_vlos*np.dot(halo_pos_vect,gal_pos_vect)
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
#            xbeta,abeta = np.loadtxt('/n/Christoq1/gifford_w/flux_code/data/average_betaprofile.tab',dtype='float',usecols=(0,1),unpack=True)
#            fit = np.polyfit((xbeta*rand_r200)[xbeta<4],abeta[xbeta<4],6)
            

            ###################################
            #End of line of sight calculations#
            ###################################
            
            #limit our particles based on set limits in projected radius, los velocity, and richness (aka: what we would observe)
            rvalues_orig,vvalues_orig,magvalues = C.set_sample(r,v,gal_mags,rlimit,vlimit,H0)
            rvalues_orig,vvalues_orig,magvalues = C.Limit_richness(rvalues_orig,vvalues_orig,magvalues,rlimit,N=rich_lim)

            #Use radial selection template (Not using now)
            #rmems,vmems = C.Selection(rmems,vmems,magmems.T[2],rand_r200)
            rnums[j] = rvalues_orig.size
            
            #create loop for bootstrap errors
            if jacknife is not None: lloop = rich_lim
            else: lloop = 1
            M_boot = np.zeros(lloop)
            for l in range(lloop):
                #randi = np.random.randint(0,rvalues_orig.size,rvalues_orig.size)
                if jacknife is None:
                    rvalues = rvalues_orig#[randi]
                    vvalues = vvalues_orig#[randi]
                else:
                    rvalues = np.append(rvalues_orig[:l],rvalues_orig[l+1:])
                    vvalues = np.append(vvalues_orig[:l],vvalues_orig[l+1:])

                gal_vdisp = astStats.biweightScale(vvalues[rvalues<rand_r200],9.0)
                rvalues = np.append(rvalues,rvalues)
                vvalues = np.append(vvalues,-vvalues)
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
                #beta = fit[0]*x_range**6 + fit[1]*x_range**5 + fit[2]*x_range**4 + fit[3]*x_range**3 + fit[4]*x_range**2 + fit[5]*x_range + fit[6]
                beta = np.zeros(x_range.size)+0.25
                g_b = (3.0-2.0*beta)/(1-beta)
                #####################
                #End Beta Definition#
                #####################

                ####
                #Do potential calculations to add to plot
                ####
                potential_NFW = np.sqrt(2*(4*np.pi*4.5e-48*(Halo_scale[i])*((Halo_srad[i]))**2*np.log(1+x_range/(Halo_srad[i]))/(x_range/(Halo_srad[i]))))*3.08e19
                pot_x,pot_y,pot_z,pot,pot_r = np.loadtxt('/n/Christoq1/MILLENNIUM/particles/halo_'+HaloID[i]+'_phi.dat',dtype='float',usecols=(0,1,2,3,5),unpack=True)
                pot_x,pot_y,pot_z,pot,pot_r = pot_x/(1+HaloZ),pot_y/(1+HaloZ),pot_z/(1+HaloZ),pot*(1+HaloZ),pot_r/(1+HaloZ)
                pot_vesc = np.sqrt(2.0*(6.67e-11*(3.24e-23)**3.0/(5.027e-31)*pot*8.6e8)*((3.0856e19)**2))
                (n,bins) = np.histogram(pot_r,200)
                pot_midbins = (bins[1:]+bins[:-1])/2.0
                here = np.digitize(pot_r,bins)
                avg_array = np.zeros(n.size)
                max_array = np.zeros(n.size)
                for uu in range(n.size):
                    try:
                        avg_array[uu] = np.average(pot_vesc[np.where(here==uu+1)])
                        max_array[uu] = np.max(pot_vesc[np.where(here==uu+1)])
                    except:
                        avg_array[uu] = 0.0
                        max_array[uu] = 0.0

                ############################
                #Caustic Surface Estimation#
                ############################
                maxv = rand_r200*H0*np.sqrt(200)+500
                Anew,threesig,dens_norm,e_dens_norm,srad,e_srad = C.level_search2(rvalues,vvalues,rvalues,vvalues,rvalues,x_range,y_range,img_tot,img_inf_tot,H0*q,HaloR200[i],rlimit,maxv,beta,Halo_srad[i],Halo_esrad[i],use_vdisp=gal_vdisp,bin=i+1)
                AnewD,threesigD,dens_normD,e_dens_normD,sradD,e_sradD = C.level_search(rvalues,vvalues,rvalues,vvalues,rvalues,x_range,y_range,img_tot,H0*q,HaloR200[i],rlimit,maxv,beta,Halo_srad[i],Halo_esrad[i],use_vdisp=gal_vdisp,bin=i+1) #use for limit richness
                vdispersion[j] = threesig/3.5
                ###########################
                #End of Caustic Estimation#
                ###########################

                M200_int = 4*np.pi*dens_norm*(srad)**3*(np.log(1+rand_r200/srad)-rand_r200/srad/(1+rand_r200/srad))
                M200_intD = 4*np.pi*dens_normD*(sradD)**3*(np.log(1+rand_r200/sradD)-rand_r200/sradD/(1+rand_r200/sradD))
                Gus_mass = 1e15*(gal_vdisp/1082.9)**(1/0.3361)
            
                #Calculate the mass profile using the caustic equation. The beta profile is optional.
                massprofile,integrand = C.masscalc(x_range,np.abs(C.Ar_finalD),rand_r200,HaloM200[i],vdispersion[j],beta=beta,conc=rand_r200/Halo_srad[i])
                massprofile3,integrand3 = C.masscalc(x_range,np.abs(C.Ar_finalD),rand_r200,HaloM200[i],vdispersion[j],beta=beta,conc=rand_r200/Halo_srad[i],dstyle=True)
                massprofile2,integrand2 = C.masscalc(x_range,np.abs(C.Ar_final),rand_r200,HaloM200[i],vdispersion[j],beta=beta,conc=rand_r200/Halo_srad[i])

                #Identify the mass within a given radius.
                lookwithin = rand_r200  #how many Mpc to print mass within
                M200_intD2 = massprofile[np.where(x_range[x_range >=0] < lookwithin)[0][-1]]
                M200_intD3 = massprofile3[np.where(x_range[x_range >=0] < lookwithin)[0][-1]]
                M200_int2 = massprofile2[np.where(x_range[x_range >=0] < lookwithin)[0][-1]]
                line_mass[j] = M200_int
                line_vdisp[j] = gal_vdisp
                print '   M_phi mass = %e' % (line_mass[j]), 'Halo mass = %e' %(HaloM200[i]),'within %.2f'%(rand_r200)
                print '   C fbe mass = %e' % (M200_intD3), 'Halo mass = %e' %(HaloM200[i]),'within %.2f'%(rand_r200)
                print '   Caus. mass = %e' % (M200_intD2), 'Halo mass = %e' %(HaloM200[i]),'within %.2f'%(rand_r200)
                print '     Gus mass = %e' % (Gus_mass), 'Halo mass = %e' %(HaloM200[i]),'within %.2f'%(rand_r200)
                print '     looping ', j+1 #print which line of sight loop we are currently on
                M_boot[l] = M200_intD2
            surf_err[j] = np.std(np.log(M_boot[np.isfinite(M_boot)]/HaloM200[i]))
            
            ####################
            #Print Mass Profile#
            ####################
#            fi = open('/nfs/christoq_ls/giffordw/Millenium/files/mprofsbig/'+str(HaloID[i])+'data/'+str(HaloID[i])+'.'+str(j)+'.'+str(rich_lim)+'mprofgals_vdisp.'+semianalytic+'.tab','w')
#            for ii in range(x_range.size):
#                fi.write(str(x_range[ii])+'\t')
#                fi.write(str(massprofile[ii])+'\t')
#                fi.write(str(mprof_true[ii])+'\t')
#                fi.write(str(C.Ar_final[ii]*np.sqrt(g_b))+'\t')
#                fi.write(str(C.Ar_finalD[ii]*np.sqrt(g_b))+'\t')
#                fi.write('\n')
#            fi.close()
            
#            li.write(str(M200_int)+'\t')
#            li.write(str(M200_intD2)+'\t')
#            li.write(str(M200_intD3)+'\t')
#            li.write(str(gal_vdisp)+'\t')
#            li.write(str(dens_norm)+'\t')
#            li.write(str(e_dens_norm)+'\t')
#            li.write(str(srad)+'\t')
#            li.write(str(e_srad)+'\t')
#            li.write('\n')
            

#        li.close()
        Mphi_hat_NFW = M200_int
        Mphi_hat_orig = M200_int2
        caustic_mass_fbeta = M200_intD3
        caustic_mass_orig = M200_intD2
        chi_sq_fit = np.sum((Anew[1:]-C.Ar_final[1:])**2/C.Ar_final[1:])
        #Plot the typical caustic diagram for this cluster/los
#        C.plotcluster(rvalues,vvalues,rvalues,vvalues,x_range,y_range,img_tot,Anew,i+1,rich_lim,maxv,img_grad_tot,img_inf_tot,beta,pot_midbins,max_array,potential_NFW,'proj')
        print 'DONE WITH CLUSTER ',i+1
        #output the important data to file. Each halo has a separate file and will need to be combined later.
#        f = open('/nfs/christoq_ls/giffordw/Millenium/files/mprofsbig/'+str(HaloID[i])+'data/'+str(HaloID[i])+'.'+str(rich_lim)+'.nideal_data.'+semianalytic+'.tab','w')
#        f.write(str(HaloM200[i])+'\t')
#        f.write(str(rand_r200)+'\t')
#        f.write(str(rvalues.size)+'\t')
#        f.write(str(Mphi_hat_NFW)+'\t')
#        f.write(str(Mphi_hat_orig)+'\t')
#        f.write(str(caustic_mass_orig)+'\t')
#        f.write(str(caustic_mass_fbeta)+'\t')
#        f.write(str(Gus_mass)+'\t')
#        f.write(str(gal_vdisp)+'\t')
#        f.write(str(dens_norm)+'\t')
#        f.write(str(e_dens_norm)+'\t')
#        f.write(str(srad)+'\t')
#        f.write(str(e_srad)+'\t')
#        f.write(str(chi_sq_fit)+'\t')
#        f.write('\n')
#        f.close()
#        if jacknife is not None:
            #ONLY PRINT WHEN DOING BOOTSTRAP ANALYSIS!!!!
#            f = open('/nfs/christoq_ls/giffordw/Millenium/files/mprofsbig/'+str(HaloID[i])+'data/'+str(HaloID[i])+'.'+str(rich_lim)+'.nideal_surferror.'+semianalytic+'.tab','w')
#            f.write(str(HaloID[i])+'\t')
#            f.write(str(np.average(surf_err))+'\t')
#            f.write('\n')
#            f.close()
        
