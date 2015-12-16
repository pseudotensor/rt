#!/bin/python

#####################################################
# 1D diagnostics for ASTRORAY
# Roman Gold, started 2015-02-01
#####################################################

import sys,os,string

import matplotlib 
#matplotlib.use("wx")
try:
    __IPYTHON__ ## ARE WE RUNNING FROM IPYTHON?
except:
    matplotlib.use("Agg") # produces png without X

import pylab,string,scipy,commands
from pylab import *
from scipy import *
from scipy.interpolate import *
from scipy.constants import *
from scipy import fftpack
import matplotlib.ticker as ticker

HOME=commands.getoutput("echo $HOME")
# RT_DIR="/rt/"
RT_DIR="/codes/rt-git/"

# assumes obs.txt in same dir (provided by Andrew Chael see [eht_python_for_roman.zip])
try:
    TELESCOPES=['ALMA','SMA','SMT','LMT','PV','PDB','SPT']

    # get all combinations
    [T1+"-"+T2 for T1 in TELESCOPES for T2 in TELESCOPES if T2!=T1]

    EHT_config_file="obs.txt"
    # EHT_config_file="obs-SMT-SMA.txt"
    eht_obs_uv = loadtxt(HOME+RT_DIR+EHT_config_file,usecols=[0,4,5],comments='#')
    fd=open(EHT_Config_file)
    EHT_config = fd.readlines()[10:] # 2nd and 3rd col: telescope sites
    fd.close()

    ALMA_BASELINES=[LINE for LINE in EHT_config if "ALMA" in LINE.split()[1:3]]
    ALMA_tuv = [[LINE.split()[0],LINE.split()[4],LINE.split()[5]] for LINE in EHT_config if "ALMA" in LINE.split()[1:3]]

    # ADD SMA,SMT,LMT,PV,PDB,SPT

except:
    pass
# scatter(obs[:,1],obs[:,2])


## USER SPECS ##
FILE_EXT="png"
mbreve="no"
dEVPA="no"
POLARIZATION_CAP = "NO" # "YES" # ARTIFICALLY CAP POLARIZATION?
m_LP_cap = 0.7
m_CP_cap = 1.0
m_LP_floor = 0.5
m_CP_floor = 0.0

FILES_2D = [FILE for FILE in sys.argv[1:] if "shotimag" in FILE]
FILES_1D = [FILE for FILE in sys.argv[1:] if "polires" in FILE or "bestfit" in FILE or "quick" in FILE or "ava" in FILE]

PLOT_SED="no"
PLOT_CORRELATED_FLUX="yes"
PLOT_I_vs_mbreve="yes"
PLOT_I_vs_vbreve="no"

if size(FILES_2D)==0:
    PLOT_CORRELATED_FLUX="no"
    PLOT_I_vs_mbreve="no"
    PLOT_I_vs_vbreve="no"

print "Found ",size(FILES_1D),"1d files and ",size(FILES_2D),"2d files\n"
print "MODES: PLOT_SED,PLOT_CORRELATED_FLUX,PLOT_I_vs_mbreve,POLARIZATION_CAP,mbreve,dEVPA",PLOT_SED,PLOT_CORRELATED_FLUX,PLOT_I_vs_mbreve,POLARIZATION_CAP,mbreve,dEVPA,"\n"

col_SED=[1,2]
if len(FILES_1D)>0:
    PLOT_SED="yes"
    SED_SCS_filename=FILES_1D[0]
    if "bestfit" in FILES_1D[0]:
        bestfit="yes"
        col_SED=[0,1,2,3,4]
    elif "polires" in FILES_1D[0]:
        polires="yes"
        col_SED=[1,2,3,4,5] # 1:nu, 2: F, 3: LP, 4: EVPA, 5: CP
    elif "quick" in FILES_1D[0]:
        polires="yes"
        col_SED=[0,1,2,3,4] # 1:nu, 2: F, 3: LP, 4: EVPA, 5: CP
    
# try:
#     SED_filenames=[element for element in sys.argv if "polires" in element]
#     # SED_SCS_filename=sys.argv[ sys.argv.index() ]"/home/rgold/rt/thickdisk7/Te35-SCS/avg/poliresa93.75th140fn6130hi.dat" # WIP
#     SED_SCS_filename=SED_filenames[0]
#     PLOT_SED="YES" # if no file found above -> IndexError

# except IndexError:
#     PLOT_SED="NO"

rc('font',size=20)

# SCATTERING could have an effect on these diagnostics...
SCATTERING = "ON" # "OFF" "ON"
if SCATTERING == "ON":
    try:
        from scipy import ndimage
    except:
        pass

angle_unit="arcsec" # "rad"
# This is used for enhancing smoother uv data at large scales: Warning not good for mtilde because the xy-data don't seem to decay to zero for large xy. Zeropadding then introducing strong artefacts.
zeropadding_factor=4

pc = scipy.constants.parsec # SI
G = scipy.constants.G # SI
c = scipy.constants.c # SI

#t_rg =  # time in secs assuming Sgr A* 
Msun = 2e30 # SI
M = 4.3e6 * Msun # SAG A*
rg = G*M/c**2
d_SagA = 8.3e3*pc
rad2microarcsec = 360/(2*pi)*3600*1e6

observing_frequency = 230. # in Ghz
print "HARDCODED frequency=",observing_frequency,"Ghz"

image_size_astroray = 8.+(600./observing_frequency)**1.5 # @230Ghz see sftab array in [ASTRORAY_main.cpp]
#image_size_astroray = 12.2 # @230Ghz see sftab array in [ASTRORAY_main.cpp]
image_size = image_size_astroray * (2*rg)/d_SagA * rad2microarcsec
image_size_rad = image_size_astroray * (2*rg)/d_SagA
shadow_schwarzschild = sqrt(27)*2*rg # diameter
shadow_maximally_spinning = 9./2.*2*rg # diameter
#######################################################


t           = [] # to be filled with simulation time
mbreve_vs_t = [] # to be filled with single point in (u,v) data
mbreve_opposite_vs_t = []
dEVPA_vs_t  = [] # to be filled with single/two opposite points in (u,v) data

# SORT IT ACCORDING TO TIME STAMP IN FILENAME
try:
    FILES_1D.sort(key=lambda elem: float(elem.split("fn")[1][:-3].split('_')[0].split('hi')[0]))
    FILES_2D.sort(key=lambda elem: float(elem.split("fn")[1][:-3].split('_')[0]))
except:
    pass

for snapshot in FILES_2D:

    print snapshot
    filename=snapshot
    fp = open(filename,"rb")
    header = fromfile(fp,count=20)
    nxy=int(header[2])+1
    data = fromfile(fp,dtype=float64).reshape(nxy,nxy,5) 
    fp.close()

    pixeldim = image_size/nxy # Specify the linear size of a pixel, in \[Mu]as

    X = pixeldim*arange(-round(nxy/2),round(nxy/2)+1)
    #?X = pixeldim*arange(-round(nxy/2),round(nxy/2))
    Y = X[:]

    freq_unit=1e-9 # uv plane scale
    uvspacing = image_size_rad/nxy

    # 150microarcsec = 1.4 in u-v plane 
    # 15 microarcsec = 14  in u-v plane
    # for plot labels given the
    # nominal scaling of 10^9\lambda for the baseline: shadow size .ie. diameter
    uv_schwarzschild = 150.*1.4/(shadow_schwarzschild /d_SagA * rad2microarcsec)
    uv_maximally_spinning = 150.*1.4/(shadow_maximally_spinning /d_SagA * rad2microarcsec)

    if string.lower(SCATTERING)=="on" or string.lower(SCATTERING)=="yes":
        FWHM = 1e3*((c/observing_frequency*freq_unit)/1e-2)**2. # valid for general frequency
        # sigma_in_microarcsec = 7.4 # valid for f=230Ghz
        sigma_in_microarcsec = FWHM/2.3

        # TESTING
        # sigma_in_microarcsec *= 5. 
        # sigma_in_microarcsec *= 2. 
        sigma_in_microarcsec *= 0.5 # Scattering can be partially undone, so... REF: http://adsabs.harvard.edu/abs/2014ApJ...795..134F


        #sigma_in_microarcsec = 1. # mild scattering to see jet but maybe kill off spurious pole features
        sigma = sigma_in_microarcsec / (image_size/nxy) # 10. # pixel units for ndimage.gaussian_filter
        for CHANNEL in range(shape(data)[-1]):
            data[:,:,CHANNEL] = ndimage.gaussian_filter(data[:,:,CHANNEL],sigma=sigma) ## sigma is in pixel units
            
    I_xy = data[:,:,0]
    Q_xy = data[:,:,1]
    U_xy = data[:,:,2]
    V_xy = data[:,:,3]

    def m_LP(I,Q,U):
        '''Compute fractional linear polarization'''
        return sqrt(Q**2+U**2)/I
    def m_CP(I,V): 
        '''Compute fractional circular polarization'''
        return abs(V)/I

############ MESS WITH POLARIZATION DATA ##############
    # if string.lower(UV_mbreve_CAP)=="yes":
    #     UV_mask  = sqrt(UV[0]**2 + UV[1]**2) > 4.
    #     UV_mask *= sqrt(UV[0]**2 + UV[1]**2) < 2.5
    #     UV_mask *= mbreve_uv > m_LP_cap
    #     UV_mask *= mbreve_uv < m_LP_floor

    if string.lower(POLARIZATION_CAP)=="yes":
        # m_LP_mask = (m_LP(I_xy,Q_xy,U_xy) < m_LP_cap)  * (m_LP(I_xy,Q_xy,U_xy) > m_LP_floor)
        m_LP_mask = (m_LP(I_xy,Q_xy,U_xy) > m_LP_cap) * (m_LP(I_xy,Q_xy,U_xy) < m_LP_floor)
        m_CP_mask = m_CP(I_xy,V_xy) > m_CP_cap

        m_LP_original = m_LP(I_xy[m_LP_mask],Q_xy[m_LP_mask],U_xy[m_LP_mask])
        m_CP_original = m_CP(I_xy[m_CP_mask],V_xy[m_CP_mask])
        Q_xy[m_LP_mask] *= m_LP_cap / m_LP_original
        U_xy[m_LP_mask] *= m_LP_cap / m_LP_original
        V_xy[m_CP_mask] *= m_CP_cap / m_CP_original
#######################################################

    I_uv = fftpack.fftshift(fftpack.fft2(I_xy,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
    Q_uv = fftpack.fftshift(fftpack.fft2(Q_xy,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
    U_uv = fftpack.fftshift(fftpack.fft2(U_xy,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
    V_uv = fftpack.fftshift(fftpack.fft2(V_xy,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))

    u = unique(fftfreq(shape(I_uv)[0],d=uvspacing)*freq_unit)
    v = unique(fftfreq(shape(I_uv)[1],d=uvspacing)*freq_unit)
    UV = meshgrid(u,v)

    mbreve_uv = abs((Q_uv+1j*U_uv)/I_uv)
    vbreve_uv = abs(V_uv)/abs(I_uv)
    EVPA_uv = angle((Q_uv+1j*U_uv)/I_uv) * 90./pi 


########### MESS WITH POLARIZATION UV DATA ###########
    UV_mbreve_CAP = "NO"
    if string.lower(UV_mbreve_CAP)=="yes":
        UV_mbreve_mask  = sqrt(UV[0]**2 + UV[1]**2) > 4.
        V_uv[UV_mbreve_mask] = 0
        UV_mbreve_mask  = sqrt(UV[0]**2 + UV[1]**2) < 2.5
        V_uv[UV_mbreve_mask] = 0
        UV_mbreve_mask  = mbreve_uv > m_LP_cap
        V_uv[UV_mbreve_mask] = 0
        UV_mbreve_mask  = mbreve_uv < m_LP_floor
        V_uv[UV_mbreve_mask] = 0

        V_uv[UV_mbreve_mask] = 0
######################################################

    # time = commands.getoutput("head -1 fieldline.*.bin").split()[0]
    #dt = 4. ## thickdisk7
    #t_ref = 6100 ## thickdisk7
    dt_GRMHD = 5. ## a0mad
    t_ref = 5500 # 2000
    print "[HARDWIRE-WARNING]: dt,t_ref=",dt_GRMHD,t_ref
    t += [(float(filename.split("fn")[1].split('_')[0].split('case')[0]) - t_ref)*dt_GRMHD * (G*M/c**3) /60./60.] # t in [hours]

    # Given time in hr    NOW PICK UV point (along EHT uv tracks
    try:
        uv_time_idx = pylab.find(eht_obs_uv[:,0]>=t[-1])[0]
    except:
        pass

    #RG:WIP DISTINGUISH DIFFERENT BASELINES

    u_probe=eht_obs_uv[uv_time_idx,1]
    v_probe=eht_obs_uv[uv_time_idx,2]

    ###################################################################
    # PICK A POINT
    # u_probe,v_probe=3.,3.

    # SHOULD REALLY JUST INTERPOLATE... E.G.
    # interp2d(u,v,mbreve_uv)(3,3)

    # lower_bound=2.5;upper_bound=3.
    try: # current EHT uv range 3.2 based on Fig.1 in ordered fields draft
        # lower_bound=0.5;upper_bound=1.0
        lower_bound=3.0;upper_bound=3.5
        # sin(pi/4) & cos(pi/4) so that we probe radius r=sqrt(u^2+v^2) where u=v
        u[(u>=cos(pi/4.)*lower_bound)*(u<=cos(pi/4.)*upper_bound)][0]
        v[(v>=sin(pi/4.)*lower_bound)*(v<=sin(pi/4.)*upper_bound)][0]
    except IndexError:
        print "No point in interval (",lower_bound,"<u<",upper_bound,"): Adjust bounds!"
        exit()

    u_probe_index = list(u).index(u[(u>lower_bound)*(u<upper_bound)][0])
    v_probe_index = list(v).index(v[(v>lower_bound)*(v<upper_bound)][0])
    u_probe_opposite_index = list(u).index(u[(u<-lower_bound)*(u>-upper_bound)][0])
    v_probe_opposite_index = list(v).index(v[(v<-lower_bound)*(v>-upper_bound)][0])
    ###################################################################

    # NEW WAY ALONG EHT TRACKS
    mbreve_vs_t += [interp2d(u,v,mbreve_uv)(u_probe,v_probe)]
    mbreve_opposite_vs_t += [interp2d(u,v,mbreve_uv)(-u_probe,-v_probe)]
    dEVPA_vs_t   += [interp2d(u,v,EVPA_uv)(u_probe,v_probe)-interp2d(u,v,EVPA_uv)(-u_probe,-v_probe)]

    # time = commands.getoutput("head -1 fieldline.*.bin").split()[0]
    dt_GRMHD = 4. ## thickdisk7
    #t_ref = 6100 ## thickdisk7
    # dt_GRMHD = 5. ## a0mad
    t_ref = 2000
    print "[HARDWIRE-WARNING]: dt,t_ref=",dt_GRMHD,t_ref
    t += [(float(filename.split("fn")[1].split('_')[0].split('case')[0]) - t_ref)*dt_GRMHD * (G*M/c**3) /60./60.]

    # OLD HARDCODED WAY
    # mbreve_vs_t += [mbreve_uv[u_probe_index,v_probe_index]]
    # mbreve_opposite_vs_t += [mbreve_uv[u_probe_opposite_index,v_probe_opposite_index]]
    # dEVPA_vs_t   += [EVPA_uv[u_probe_index,v_probe_index]-EVPA_uv[u_probe_opposite_index,v_probe_opposite_index]]
    # # dEVPA_vs_t   += [-EVPA_uv[u_probe_index,v_probe_index]+EVPA_uv[u_probe_opposite_index,v_probe_opposite_index]]

for jump_removal_iteration in range(5):
    try:
        jump_down_in_EVPA = dEVPA_vs_t.index(array(dEVPA_vs_t)[diff(dEVPA_vs_t)<-100][0])
        dEVPA_vs_t[jump_down_in_EVPA+1:] = array(dEVPA_vs_t)[jump_down_in_EVPA+1:] + 180.
    except:
        pass
    try:
        jump_up_in_EVPA = dEVPA_vs_t.index(array(dEVPA_vs_t)[diff(dEVPA_vs_t)>100][0])
        dEVPA_vs_t[jump_up_in_EVPA+1:] = array(dEVPA_vs_t)[jump_up_in_EVPA+1:] - 180.
    except:
        pass



if size(FILES_2D)>0:

  if mbreve=="yes":
      ##########################
      figure(0,figsize=(12,6)) #

      labelstring_mbreve=[
          r"$uv="+str(round(u[u_probe_index],1))+"G\lambda$",
          r"$uv="+str(round(u[u_probe_opposite_index],1))+"G\lambda$"
          # r"$(u="+str(round(u[u_probe_index],1))+",v="+str(round(v[v_probe_index],1))+")$",
          # r"$(u="+str(round(u[u_probe_opposite_index],1))+",v="+str(round(v[v_probe_opposite_index],1))+")$"
          ]
      plot(t,mbreve_vs_t,"bo-",label=labelstring_mbreve[0])
      plot(t,mbreve_opposite_vs_t,"r+-",label=labelstring_mbreve[1])
      legend(loc="upper right",labelspacing=0.1) # ,fontsize=15)
      title(os.getcwd().split('/')[-1])
      xlabel(r"$t/h$")
      ylabel(r"$\breve{m}$")
      # ylabel(r"$\breve{m}(\|u\|="+str(round(v[v_probe_index],1))+",\|v\|="+str(round(v[v_probe_index],1))+")$")
      axis((amin(t),amax(t),0,1.1))
      #axis((t[0],t[-1],0,1.1))
      tight_layout()
      savefig("mbreve-vs-t.png")
      savefig("mbreve-vs-t.pdf")


  if dEVPA=="yes":
      ########################
      figure(1,figsize=(12,6)) #

      # WIP clean 180deg jumps using diff()
      plot(t,dEVPA_vs_t,"g^")
      xlabel(r"$t/h$")
      ylabel(r"$\delta EVPA(\|uv\|="+str(round(v[v_probe_index],1))+"G\lambda)$") # +",\|v\|="+str(round(v[v_probe_index],1))+")$")
      #xlim(t[0],t[-1])
      tight_layout()
      savefig("dEVPA-vs-t.png")
      savefig("dEVPA-vs-t.pdf")

  if False:
      figure(2)
      # FIXME: still need to filter |uv|= 2.5-4
      # CHECK THIS
      hist( (abs(V_uv)/abs(I_uv))[(mbreve_uv < m_LP_cap) * (mbreve_uv > m_LP_floor)].flatten() ,50,normed=True)
      # CHECK THIS
      CP_mask = (mbreve_uv < m_LP_cap) * (mbreve_uv > m_LP_floor) * ( sqrt(UV[0]**2 + UV[1]**2) < 4.) * ( sqrt(UV[0]**2 + UV[1]**2) > 2.5 )
      hist( (abs(V_uv)/abs(I_uv))[(mbreve_uv < m_LP_cap) * (mbreve_uv > m_LP_floor)].flatten() ,50,normed=True)
      # CHECK THIS
      hist( (abs(V_uv)/abs(I_uv))[CP_mask].flatten() ,19)

      figure(3)
      # V_uv[CP_mask] = 0
      # V_uv[UV_mbreve_mask] = 0
      pcolormesh(UV[0],UV[1],abs(V_uv)/abs(I_uv),cmap=cm.bone,vmax=None)
      axis((-5,5,-5,5))
      colorbar(pad=0)
      title(r"$\rm \|\tilde{V}\|/\|\tilde{I}\|$ where $0.5<\breve{m}<0.7$") # replace with parameters!
      ylabel(r"$v/G\lambda$")
      xlabel(r"$u/G\lambda$")
      tight_layout()



#######################################
################ SED ##################

if string.lower(PLOT_SED)=="yes":

    figure(4)
    #SED_SCS_filename="/home/rgold/rt/thickdisk7/Te35-SCS/avg/poliresa93.75th140fn6130hi.dat" # WIP
    SED_SCS=loadtxt(SED_SCS_filename,usecols=col_SED)
    try:
        SED_SCSjet_filename="/home/rgold/rt/thickdisk7/Te10/poliresa93.75th140fn6130hi.dat" # WIP
        SED_SCSjet=loadtxt(SED_SCSjet_filename,usecols=col_SED)
        #SED_default_filename="../../Te-default/poliresa93.75th140fn6130hi.dat" # WIP
        SED_default_filename="/home/rgold/rt/thickdisk7/Te-default/poliresa93.75th130fn6130hi.dat" # WIP
        SED_default=loadtxt(SED_default_filename,usecols=col_SED)

        #plot(SED_SCSjet[:,0],SED_SCSjet[:,1],'mx-',label=r"$\rm T_{e,jet}=10m_ec^2,SCS+jet$")
        #plot(SED_default[:,0],SED_default[:,1],'cD-',label=r"$\rm default$")

    except IOError:
        pass
        
    # PLOT SED & COMPARE TO OBSERVATIONS
    nu_obs = array([8.45, 14.90, 22.50, 43.00, 87.73, 102., 145., 230.86, 349., 674., 857., 1500., 3000., 5000.])
    SED_errors = array([0.031, 0.012, 0.015, 0.026, 0.080, 0.1517, 0.2644, 0.1414, 0.1205, 0.3508, 0.2404, 0. , 0., 0.])


#polarized spectrum of Sgr A*, each array element is (frequency, Fnu, LP, EVP\A, CP)
#Fnu:  flux at frequency nu
#LP:   linear polarization fraction at frequency nu
#EVPA: Electric vector position angle
#CP:   circular polarization fraction at frequency nu
    # tofit = [[8.450, 0.683, 0., 0., -0.2500], [14.90, 0.871, 0., 0., -0.6200], [22.50, 0.979, 0.1900, 131.0, 0.], [43.00, 1.135, 0.5500, 94.25, 0.], [87.73, 1.841, 1.420, -4., 0.], [102.0, 1.908, 0., 0., 0.], [145.0, 2.275, 0., 0., 0.], [230.9, 2.637, 7.398, 111.5, -1.200], [349.0, 3.181, 6.499, 146.9, -1.500], [674.0, 3.286, 0., 0., 0.], [857.0, 2.867, 0., 0., 0.], [1500., 1., 0., 0., 0.], [3000., 1., 0., 0., 0.], [5000., 1., 0., 0., 0.]];
    tofit = [[8.450, 0.683, None, None, -0.2500], [14.90, 0.871, None, None, -0.6200], [22.50, 0.979, 0.1900, 131.0, None], [43.00, 1.135, 0.5500, 94.25, None], [87.73, 1.841, 1.420, -4., None], [102.0, 1.908, None, None, None], [145.0, 2.275, None, None, None], [230.9, 2.637, 7.398, 111.5, -1.200], [349.0, 3.181, 6.499, 146.9, -1.500], [674.0, 3.286, None, None, None], [857.0, 2.867, None, None, None], [1500., 1., None, None, None], [3000., 1., None, None, None], [5000., 1., None, None, None]];
    tofit=array(tofit)

    SgrA_SED_Observed=[[],[]] # [nu_Ghz,F_nu] [WIP]
    SgrA_LP_Observed=[[],[]] # [nu_Ghz,LP] [WIP]
    SgrA_CP_Observed=[[],[]] # [nu_Ghz,CP] [WIP]
    SgrA_EVPA_Observed=[[],[]] # [nu_Ghz,EVPA] [WIP]
    nu=linspace(0,1000,100)

    def SgrA_SED_FIT(nu):
        return 0.248*nu**0.45 * exp(-(nu/1100.)**2)


    def ChiSq(F,tofit=tofit,SED_errors=SED_errors):
        Chisq = sum((F[4:-4]-tofit[4:-4,1])**2/SED_errors[4:-4]**2)
        dof   = 7 ## 10-3: fit for i,rho and T_e
        return Chisq/dof


    #plot(SED_SCS[:,0],SED_SCS[:,1],'bo-',label=r"$\rm T_{e,jet}=35m_ec^2,SCS$")
    errorbar(nu_obs,SgrA_SED_FIT(nu_obs),yerr=SED_errors,fmt='ks',label="observed")
    #plot(nu[:-3],SgrA_SED_FIT(nu[:-3]),'k--',alpha=0.5,lw=4,label="FIT") 

    # plot_styles = ["rd-","b+"]
    # taken from orbits.py
    # style = Line2D.markers.keys()[index]+Line2D.lineStyles.keys()[index]
    # plot_styles = ["rd-","b+"]
    for FILE in FILES_1D:

        # be cyclic in finite marker arrays
        size_markers = size(Line2D.markers.keys())
        size_lineStyles = size(Line2D.lineStyles.keys())

        # taken from orbits.py
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~vv [issue with "0"]
        plot_style = str(Line2D.markers.keys()[(FILES_1D.index(FILE)+1)%size_markers])+str(Line2D.lineStyles.keys()[(FILES_1D.index(FILE)+3)%size_lineStyles])
        # plot_style = str(Line2D.markers.keys()[(FILES_1D.index(FILE)+1)%size_markers])+str(Line2D.lineStyles.keys()[(FILES_1D.index(FILE)+3)%size_lineStyles])

        # plot(SED[:,0],SED[:,1],plot_style,label=FILE)# r"$\rm T_{e,jet}=35m_ec^2,SCS$")
        if "quick" in FILE:
            SED=loadtxt(FILE,usecols=[0,1])
            plot(SED[:,0],SED[:,1],'b-',linewidth=2)
        else:
            SED=loadtxt(FILE,usecols=col_SED)
            plot(SED[:,0],SED[:,1],'c.-',alpha=0.5,label=FILE)# r"$\rm T_{e,jet}=35m_ec^2,SCS$")

        if FILES_1D==FILE:
            legend(loc="lower right",labelspacing=0.2) # ,fontsize=15)

    xlabel(r"$\nu/Ghz$")
    ylabel(r"$F_\nu/Jy$")
    xlim(0,1000) # ,1600)
    tight_layout()
    savefig("SED.png")
    savefig("SED.pdf")

if string.lower(PLOT_SED)=="yes":
    FILE=[FILE for FILE in FILES_1D if "bestfit" in FILE or "ava" in FILE or "quick" in FILE or "poli" in FILE][0]
    SED=loadtxt(FILE)

    figure(5)
    bestfit_labels=[r"$F_\nu$",r"$<LP>$",r"$<CP>$",r"$<EVPA>$"]
    yticks_obs=[arange(5),arange(0,20,5),arange(-2,5),arange(0,200,45)]
    # bestfit_labels=[r"$F_\nu/Jy$",r"$<LP>/\%$",r"$<EVPA>/{}^\circ$",r"$<CP>/\%$"]
    # yticks_obs=[arange(5),arange(0,20,5),arange(0,200,45),arange(-2,5)]
    for panel in range(1,5):
        subplot(220+panel)
        plot(SED[:,col_SED[0]],SED[:,col_SED[0]+panel],label="model")
        if panel==1:
            plot(nu,SgrA_SED_FIT(nu),'k--',alpha=0.5,lw=4,label="FIT") 
            #errorbar(nu_obs,SgrA_SED_FIT(nu_obs),yerr=SED_errors,fmt='ks',label="observed")
            errorbar(tofit[:,0],tofit[:,panel],yerr=SED_errors,fmt='ks',label="observed")
        #elif panel==4:
            # not enough space
            # legend(loc="lower right") # ,fontsize=12)
        else:
            plot(tofit[:,0],tofit[:,[None,1,2,4,3][panel]],"ks",label="observations")
        #errorbar(tofit[:,0],tofit[:,panel],"rs",yerr=SED_errors,fmt='ks',label="observed")
        gca().set(xlabel=r"$\nu$",ylabel=bestfit_labels[panel-1],xticks=arange(0,1250,250),yticks=yticks_obs[panel-1])
        xlim(None,1e3)

        tight_layout()
        savefig("F_LP_CP_EVPA.png")



####################################################
## CORRELATED FLUX DENSITY / EMISSION REGION SIZE ##

# Doeleman et al 2008 : F=(2.4+/-0.5)Jy |uv|=3.5Glambda f=?Ghz FWHM=(43+14/-8) muarcsec FWHM-scatter=(37+16/-10) muarcsec errors are 3sigma
# baselines range: 0.5Glambda - 3.5Glambda
# F(3.5Glambda)=0.35Jy simultaneously with F(0)=2.4Jy\pm0.5 at 230Ghz

if PLOT_CORRELATED_FLUX=="yes":
    figure(6)

    uv_idx = shape(I_uv)[0]/2
    
    ## interpolate on circle (so that we can sum over all angles

    r_uv_data = sqrt(u**2+v**2)
    ## ph_uv_data = arctan2(u,v)
    ph_uv_data = arctan2(v,u)

    nr=51;nph=51
    r_uv  = linspace(0,10,nr)
    ph_uv = linspace(0,2.*pi,nph)

    I_uv_max = amax(abs(I_uv))
    I_uv_intp = interp2d(u,v,abs(I_uv))
    mbreve_uv_intp = interp2d(u,v,abs(mbreve_uv))
    vbreve_uv_intp = interp2d(u,v,abs(V_uv)/abs(I_uv))

    # CF = empty((nr,nph))
    mbreve_uv_rphi,I_uv_rphi,vbreve_uv_rphi = empty((nr,nph)),empty((nr,nph)),empty((nr,nph))
    for r_idx in arange(nr):
        for phi_idx in arange(nph):
            I_uv_rphi[r_idx][phi_idx]      = I_uv_intp(r_uv[r_idx]*sin(ph_uv[phi_idx]),r_uv[r_idx]*cos(ph_uv[phi_idx]))
            mbreve_uv_rphi[r_idx][phi_idx] = mbreve_uv_intp(r_uv[r_idx]*sin(ph_uv[phi_idx]),r_uv[r_idx]*cos(ph_uv[phi_idx]))
            vbreve_uv_rphi[r_idx][phi_idx]      = vbreve_uv_intp(r_uv[r_idx]*sin(ph_uv[phi_idx]),r_uv[r_idx]*cos(ph_uv[phi_idx]))
    I_uv_1d = mean(abs(I_uv_rphi),axis=1)
    I_uv_1d_min = amin(abs(I_uv_rphi),axis=1)
    I_uv_1d_max = amax(abs(I_uv_rphi),axis=1)
    mbreve_uv_1d = mean(abs(mbreve_uv_rphi),axis=1)
    mbreve_uv_1d_min = amin(abs(mbreve_uv_rphi),axis=1)
    mbreve_uv_1d_max = amax(abs(mbreve_uv_rphi),axis=1)
    vbreve_uv_1d = mean(abs(vbreve_uv_rphi),axis=1)

    plot(r_uv,I_uv_1d/I_uv_1d[0],'cs-',label=r"$I:\phi-avg$")
    fill_between(r_uv,I_uv_1d_min/I_uv_1d[0],I_uv_1d_max/I_uv_1d[0],color='cyan',alpha=0.5)

    plot(r_uv,mbreve_uv_1d,'rd-',label=r"$\breve{m}:\phi-avg$")
    fill_between(r_uv,mbreve_uv_1d_min,mbreve_uv_1d_max,color='red',alpha=0.2)

    plot(r_uv,vbreve_uv_1d,'o-',color='gray',label=r"$\breve{v}:\phi-avg$")
    # plot(v,abs(I_uv[uv_idx,:])/I_uv_max,'kx-',label=r"$I(u=0)$")
    # plot(u,abs(I_uv[:,uv_idx])/I_uv_max,'m+-',label=r"$I(v=0)$")

    I_1d_uv_obs = array([[0.6,1],[2.8,0.1],[3,0.2],[3.5,0.2]]) #RG: ...WIP... BY EYE FROM SCIENCE PLOT
    #FIXME I_1d_uv_obs = array([[0.6,1],[2.8,None],[3,None],[3.5,0.35]])
    I_uv_err = array([0.2,0.03,0.05,0.05]) #RG: ...WIP... BY EYE FROM SCIENCE PLOT
    errorbar(I_1d_uv_obs[:,0],I_1d_uv_obs[:,1],yerr=I_uv_err,fmt='ks',label="observed (day 80)")
    # ylabel(r"$\|\tilde{I}/\tilde{I}_{\rm max}\|$");
    axis((0,10,0,1.1));xlabel(r"$\|uv\|$");legend=legend(labelspacing=0.1);tight_layout()
    legend.get_frame().set_alpha(0.5)
    plt.setp(gca().get_legend().get_texts(), fontsize='18')

    if size(FILES_2D)>0:
        iter=FILES_2D[0].split("fn")[1].split("_")[0].split("case")[0]
        # t_ref=0.
        # dt_GRMHD # ?
        t=(float(iter)-t_ref)*dt_GRMHD * (G*M/c**3) /60./60.  # in hours
        title(r"$t="+str(round(t,1))+"$h")

    tight_layout()
    savefig("Fnu-vs-uv_"+string.zfill(iter,4)+"."+FILE_EXT)


if string.lower(PLOT_I_vs_mbreve)=="yes":
    figure(8)
    scatter(abs(I_uv)/I_uv_max,mbreve_uv,c="r",marker="x",alpha=0.5)
    axis((0,1,0,1.2))
    xlabel(r"$\|I(u,v)/I_{\rm max}\|$")
    ylabel(r"$\breve{m}$")
    tight_layout()
    savefig("I-vs-mbreve."+FILE_EXT)

if string.lower(PLOT_I_vs_vbreve)=="yes":
    figure(9)
    scatter(abs(I_uv)/I_uv_max,vbreve_uv,c="m",alpha=0.5)
    axis((0,1,0,1.2))
    xlabel(r"$\|I(u,v)/I_{\rm max}\|$")
    ylabel(r"$\breve{v}$")
    tight_layout()
    savefig("I-vs-vbreve."+FILE_EXT)

# See divenex's answer on
# http://stackoverflow.com/questions/6163334/binning-data-in-python-with-scipy-numpy
# import numpy as np
# from scipy.stats import binned_statistic
# data = np.random.rand(100)
# bin_means = binned_statistic(data, data, bins=10, range=(0, 1))[0]

##################################
print "============\n====DONE===="
##################################
