'''finding masses of Millenium clusters'''
import numpy as np
from flux_caustics import *
#import pyfits
import sys
import os
import astStats
from matplotlib.pyplot import *
from scipy import weave
from scipy.weave import converters

class get_data:
    def get_particles(self,ID,H0):
        "Specify the directory/filename(s) of the particles in the simulation"
        '''
        try:
            r,part_x,part_y,part_z,part_vx,part_vy,part_vz,id = np.loadtxt('minimill/halo_'+ID+'_d35_particles.csv',dtype='float',delimiter=',',usecols=(0,1,2,3,4,5,6,7),unpack=True)
        '''
        #except:
        part_x,part_y,part_z,part_vx,part_vy,part_vz,id = np.loadtxt('minimill/minidata/halo_'+ID+'_d35_particles.csv',dtype='float',delimiter=',',usecols=(0,1,2,3,4,5,6),unpack=True)
        return (part_x/(H0/100.0),part_y/(H0/100.0),part_z/(H0/100.0),part_vx,part_vy,part_vz)

    def get_clusters(self,H0):
        "Specify the directory/filename of the halos in the simulation"
        M200,clus_vdisp,r200,cluster_xpos,cluster_ypos,cluster_zpos,cluster_vx,cluster_vy,cluster_vz = np.loadtxt('minimill/halos.csv',dtype='float',delimiter=',',usecols=(4,16,7,12,13,14,17,18,19),unpack=True)
        haloid,subhaloid = np.loadtxt('minimill/halos.csv',dtype='string',delimiter=',',usecols=(1,2),unpack=True)
        halo_conc = np.loadtxt('minimill/mini_concentration.csv',delimiter=',',dtype='float',usecols=(1,),unpack=True)
        #ang_distances = np.loadtxt('data/ang_distances_lowz2.tab',dtype='float',usecols=(0,),unpack=True)
        return (haloid,M200,subhaloid,cluster_xpos/(H0/100.0),cluster_ypos/(H0/100.0),cluster_zpos/(H0/100.0),cluster_vx,cluster_vy,cluster_vz,r200,clus_vdisp,halo_conc)

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
    H0 = 72.0
    q = 10.0
    c = 300000.0

    #load halo data from file by calling the get_clusters function
    print 'GETTING CLUSTERS'
    haloid,M200,subhaloid,clus_xpos,clus_ypos,clus_zpos,clus_vx,clus_vy,clus_vz,r200,clus_vdisp,clus_conc = G.get_clusters(H0)
    halo_p = np.array([clus_xpos,clus_ypos,clus_zpos])
    halo_v = np.array([clus_vx,clus_vy,clus_vz])
    N = haloid.size
    
    #Both False: Use Dan's sigma clipping -- 
    #Use_vdisp = True: Use the table vdisp -- 
    #use_mems = True: Use known members for vdisp (outdated)
    use_mems = False
    use_vdisp = True
    
    line_num = 50 #This specifies how many lines of sight for each cluster
    rich_lim = 100 #This specifies hwo many particles will be used to estimate the caustic surface

    #All my data files are sorted by M200_crit, so here is the sorting block
    HaloID = haloid[np.argsort(M200)[::-1]][:N]
    HaloR200 = r200[np.argsort(M200)[::-1]][:N]/(H0/100.0)
    HaloM200 = M200[np.argsort(M200)[::-1]][:N]/(H0/100.0)
    HaloVD = clus_vdisp[np.argsort(M200)[::-1]][:N]
    Halo_PX = clus_xpos[np.argsort(M200)[::-1]][:N]
    Halo_PY = clus_ypos[np.argsort(M200)[::-1]][:N]
    Halo_PZ = clus_zpos[np.argsort(M200)[::-1]][:N]
    Halo_VX = clus_vx[np.argsort(M200)[::-1]][:N]
    Halo_VY = clus_vy[np.argsort(M200)[::-1]][:N]
    Halo_VZ = clus_vz[np.argsort(M200)[::-1]][:N]
    Halo_conc = clus_conc[np.argsort(M200)[::-1]][:N]

    
    try: #if the user feeds a number at the command line, that halo is used. It should usually be 0.
        ar = np.int(sys.argv[1])
    except IndexError: #if the user does not feed a number, the index is chosen to start with 0
        ar = 0
    vdisp_crit = np.zeros(16)
    for k in range(1): #loop over the halos
        #i = ar
        i = k
        mem_flag = 1 #This is a flag to alerting to (1) if members > 0 and (0) if not. Affects number output.
        #Get current galaxy info
        print 'GETTING GALAXIES'
        part_xpos,part_ypos,part_zpos,part_vx,part_vy,part_vz = G.get_particles(HaloID[i],H0)
        gal_p = np.array([part_xpos,part_ypos,part_zpos])
        gal_v = np.array([part_vx,part_vy,part_vz])

        #organize the current halo position and velocity
        Halo_P = np.array([Halo_PX[i],Halo_PY[i],Halo_PZ[i]]) #current halo position
        Halo_V = np.array([Halo_VX[i],Halo_VY[i],Halo_VZ[i]]) #current halo velocity

        #calculate the radial position of every particle we have loaded, as well as the total velocity.
        r_pos = np.sqrt((part_xpos-Halo_P[0])**2+(part_ypos-Halo_P[1])**2+(part_zpos-Halo_P[2])**2)
        
        #calculate vdisp of particles
        HVD = np.sqrt(astStats.biweightScale(part_vx[np.where(r_pos<=HaloR200[i])]-Halo_V[0],9.0)**2+astStats.biweightScale(part_vy[np.where(r_pos<=HaloR200[i])]-Halo_V[1],9.0)**2+astStats.biweightScale(part_vz[np.where(r_pos<=HaloR200[i])]-Halo_V[2],9.0)**2)/np.sqrt(3)

        #Define the filename(s) to output the results
        if use_mems == True:
            f = open('files/'+str(rich_lim)+'n/'+str(HaloID[i])+'.minimill_masses_'+str(rich_lim)+'n_mems.tab','w')
        elif use_vdisp == True:
            f = open('files/'+str(rich_lim)+'n/'+str(HaloID[i])+'.minimill_masses_'+str(rich_lim)+'n_vdisp_los_part.tab','w')
        else:
            f = open('files/'+str(rich_lim)+'n/'+str(HaloID[i])+'.minimill_masses_'+str(rich_lim)+'n_los_part.tab','w')
        
        line_mass = np.zeros(line_num) #empty array for different lines of sight masses
        vdispersion = np.zeros(line_num) #empty array for different lines of sight vdisp
        #line_error = np.zeros(line_num) #empty array for the errors above

        for j in range(line_num): #loop over different lines of sight
            #define r200 and limits to select our particles for the caustic estimation
            rand_r200 = HaloR200[i]
            rlimit = rand_r200*1.25
            vlimit = 3500.0

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
            r = (halo_dist**2-np.dot(halo_pos_vect*halo_dist,gal_pos_vect)**2)**0.5
            v = gal_vlos-halo_vlos*np.dot(halo_pos_vect,gal_pos_vect) #dot product is cosine of the angle between the halo center and galaxy
            ###################################
            #End of line of sight calculations#
            ###################################
            
            #limit our particles based on set limits in projected radius and los velocity (aka: what we would observe)
            rvalues = np.array(r)[np.where((r < rlimit) & (v < vlimit) & (v > -vlimit))]
            vvalues = np.array(v)[np.where((r < rlimit) & (v < vlimit) & (v > -vlimit))]
            try:
                vfix = astStats.biweightLocation(vvalues[np.where((rvalues<0.5) & (vvalues>-vlimit) & (vvalues<vlimit))],6.0)
                vvalues = vvalues - vfix
            except:
                vfix = np.average(vvalues[np.where((rvalues<0.5) & (vvalues>-3500) & (vvalues<3500))])
                vvalues = vvalues - vfix
            rvalues1 = np.array(rvalues)[np.where((rvalues < rlimit) & (vvalues < vlimit) & (vvalues > -vlimit))]
            vvalues1 = np.array(vvalues)[np.where((rvalues < rlimit) & (vvalues < vlimit) & (vvalues > -vlimit))]
            rvalues = rvalues1
            vvalues = vvalues1
			
            #limit our sample based on a richness limit
            rand_int = np.random.randint(0,rvalues.size,rich_lim)
            rvalues = rvalues[rand_int]
            vvalues = vvalues[rand_int]

            
            ###########################
            #Kernal Density Estimation#
            ###########################
            res_array = np.array([10,5,3,2,1])
            img_tot = 0
            for u in range(res_array.size):
                x_range,y_range,img = C.gaussian_kernel(rvalues,vvalues,rand_r200,normalization=H0,scale=q,res=200,adj = res_array[u],see=False)
                img_tot += img/np.max(img)
            ###########################
            #End of Density Estimation#
            ###########################

            
            ############################
            #Caustic Surface Estimation#
            ############################
            maxv = rand_r200*H0*np.sqrt(200)+500
            if use_vdisp == True:
                Anew,threesig = C.level_search(rvalues,vvalues,rvalues,vvalues,rvalues,x_range,y_range,img_tot,H0*q,rand_r200,rlimit,maxv,use_vdisp=HVD,use_mems=use_mems) #use for limit richness
            else:
                Anew,threesig = C.level_search(rvalues,vvalues,rvalues,vvalues,rvalues,x_range,y_range,img_tot,H0*q,rand_r200,rlimit,maxv,use_mems=use_mems) #use for limit richness
            vdispersion[j] = threesig/3.5
            ###########################
            #End of Caustic Estimation#
            ###########################
            
            #Use particle positions/velocities to calculate a beta profile for use in the mass calculation.
            print 'CALCULATING BETA'
            beta = C.betaprofile(part_xpos,part_ypos,part_zpos,part_vx,part_vy,part_vz,Halo_P[0],Halo_P[1],Halo_P[2],Halo_V[0],Halo_V[1],Halo_V[2],x_range,rlimit)
            
            #(NOT WORKING) Use particle postions/velocities to calculate a density profile for use in the mass calculation.
            #print 'CALCULATION DENSITY'
            #density,density_tot = C.densityprofile(part_xpos,part_ypos,part_zpos,Halo_P[0],Halo_P[1],Halo_P[2],x_range)
            #massprofile = C.masscalc(x_range,np.abs(Anew),rand_r200,HaloM200[i],vdispersion[j],density=density,density_tot=density_tot,beta=beta)
            
            #Calculate the mass profile using the caustic equation. The beta profile is optional.
            massprofile = C.masscalc(x_range,np.abs(Anew),rand_r200,HaloM200[i],vdispersion[j],beta=beta,conc=Halo_conc[i])
            
            #Identify the mass within a given radius.
            lookwithin = rand_r200
            line_mass[j] = massprofile[np.where(x_range[x_range >=0] < lookwithin)[0][-1]]
            print '     mass = %e' % (line_mass[j]), 'Halo mass = %e' % (HaloM200[i]*1e10),'within %.2f'%(rand_r200)
            print '     looping ', j+1 #print which line of sight loop we are currently on

            #print the mass profile out to be used later.
            try:
                fi = open('files/mprofs/'+str(HaloID[i])+'data/'+str(HaloID[i])+'.'+str(j)+'.'+str(rich_lim)+'mprof.tab','w')
            except IOError:
                print 'Creating directory for '+str(HaloID[i])+' under files/mprofs/'
                os.system('mkdir files/mprofs/'+str(HaloID[i])+'data')
                fi = open('files/mprofs/'+str(HaloID[i])+'data/'+str(HaloID[i])+'.'+str(j)+'.'+str(rich_lim)+'mprof.tab','w')
            mprof_true = np.zeros(x_range.size)
            for ii in range(x_range.size):
                fi.write(str(x_range[ii])+'\t')
                fi.write(str(massprofile[ii])+'\t')
                if ii != 0:
                    mprof_true[ii] = 8.6e8 * r_pos[np.where(r_pos<=x_range[ii])].size/(H0/100.0)
                    fi.write(str(mprof_true[ii])+'\t')
                else:
                    fi.write(str(0.0)+'\t')
                fi.write(str(Anew[ii])+'\t')
                fi.write('\n')
            fi.close()
            '''
            #print the beta profile out to be used later
            fi = open('files/mprofs/'+str(HaloID[i])+'data/'+str(HaloID[i])+'.beta.tab','w')
            for ii in range(x_range.size):
                fi.write(str(x_range[ii])+'\t')
                fi.write(str(beta[ii])+'\t')
                fi.write('\n')
            fi.close()
            '''
        
            #Plot the typical caustic diagram for this cluster/los
            #C.plotcluster(rvalues,vvalues,x_range,y_range,img_tot,Anew,i+1,rich_lim,maxv)
        print 'DONE WITH CLUSTER ',i+1
        
        #output the important data to file. Each halo has a separate file and needs to be combined later.
        f.write(str(HaloM200[i]*1e10)+'\t')
        f.write(str(rand_r200)+'\t')
        f.write(str(rvalues.size)+'\t')
        for h in range(line_mass.size):
            f.write(str(line_mass[h])+'\t')
        f.write('\n')
        f.write(str(HaloM200[i]*1e10)+'\t')
        f.write(str(rand_r200)+'\t')
        f.write(str(rvalues.size)+'\t')
        for h in range(line_mass.size):
            f.write(str(vdispersion[h])+'\t')
        #f.write('\n')
        f.close()
        vdisp_crit[i] = HVD
