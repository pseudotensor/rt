#!/bin/python

#####################################################
# 1D diagnostics for ASTRORAY
# Roman Gold, started 2015-02-01
#####################################################

import sys

import matplotlib 
#matplotlib.use("wx")
try:
    __IPYTHON__ ## ARE WE RUNNING FROM IPYTHON?
except:
    matplotlib.use("Agg") # produces png without X

import pylab,string,scipy,commands
from pylab import *
from scipy import *
from scipy.constants import *
from scipy import fftpack
import matplotlib.ticker as ticker

## USER SPECS ##
POLARIZATION_CAP = "NO" # "YES" # ARTIFICALLY CAP POLARIZATION?
m_LP_cap = 0.7
m_CP_cap = 1.0
m_LP_floor = 0.5
m_CP_floor = 0.0

try:
    SED_filenames=[element for element in sys.argv if "polires" in element]
    # SED_SCS_filename=sys.argv[ sys.argv.index() ]"/home/rgold/rt/thickdisk7/Te35-SCS/avg/poliresa93.75th140fn6130hi.dat" # WIP
    SED_SCS_filename=SED_filenames[0]
    SED="YES" # if no file found above -> IndexError

except IndexError:
    SED="NO"

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

FILES_2D = [FILE for FILE in sys.argv[1:] if "shotimag" in FILE]
FILES_1D = [FILE for FILE in sys.argv[1:] if "polires" in FILE]

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
        sigma_in_microarcsec = 7.4 # valid for f=230Ghz
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

    ###################################################################
    # PICK A POINT
    # u_probe,v_probe=3.,3.
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

    # time = commands.getoutput("head -1 fieldline.*.bin").split()[0]
    t += [(float(filename.split("fn")[1].split('_')[0])-6100.)*4. * (G*M/c**3) /60.]
    mbreve_vs_t += [mbreve_uv[u_probe_index,v_probe_index]]
    mbreve_opposite_vs_t += [mbreve_uv[u_probe_opposite_index,v_probe_opposite_index]]
    dEVPA_vs_t   += [EVPA_uv[u_probe_index,v_probe_index]-EVPA_uv[u_probe_opposite_index,v_probe_opposite_index]]
    # dEVPA_vs_t   += [-EVPA_uv[u_probe_index,v_probe_index]+EVPA_uv[u_probe_opposite_index,v_probe_opposite_index]]

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


  ##########################
  figure(0,figsize=(12,6)) #

  labelstring_mbreve=[
      r"$uv="+str(round(u[u_probe_index],1))+"G\lambda$",
      r"$uv="+str(round(u[u_probe_opposite_index],1))+"G\lambda$"
      # r"$(u="+str(round(u[u_probe_index],1))+",v="+str(round(v[v_probe_index],1))+")$",
      # r"$(u="+str(round(u[u_probe_opposite_index],1))+",v="+str(round(v[v_probe_opposite_index],1))+")$"
      ]
  plot(t,mbreve_vs_t,"bo",label=labelstring_mbreve[0])
  plot(t,mbreve_opposite_vs_t,"r+",label=labelstring_mbreve[1])
  legend(loc="upper right",labelspacing=0.1,fontsize=15)
  xlabel(r"$t/min$")
  ylabel(r"$\breve{m}$")
  # ylabel(r"$\breve{m}(\|u\|="+str(round(v[v_probe_index],1))+",\|v\|="+str(round(v[v_probe_index],1))+")$")
  axis((amin(t),amax(t),0,1.1))
  #axis((t[0],t[-1],0,1.1))
  tight_layout()
  savefig("mbreve-vs-t.png")


  ########################
  figure(1,figsize=(12,6)) #

  # WIP clean 180deg jumps using diff()
  plot(t,dEVPA_vs_t,"g^")
  xlabel(r"$t/min$")
  ylabel(r"$\delta EVPA(\|uv\|="+str(round(v[v_probe_index],1))+"G\lambda)$") # +",\|v\|="+str(round(v[v_probe_index],1))+")$")
  #xlim(t[0],t[-1])
  tight_layout()
  savefig("dEVPA-vs-t.png")

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

if string.lower(SED)=="yes":

    figure(4)
    
    #SED_SCS_filename="/home/rgold/rt/thickdisk7/Te35-SCS/avg/poliresa93.75th140fn6130hi.dat" # WIP
    SED_SCS=loadtxt(SED_SCS_filename,usecols=(1,2))
    try:
        SED_SCSjet_filename="/home/rgold/rt/thickdisk7/Te10/poliresa93.75th140fn6130hi.dat" # WIP
        SED_SCSjet=loadtxt(SED_SCSjet_filename,usecols=(1,2))
        #SED_default_filename="../../Te-default/poliresa93.75th140fn6130hi.dat" # WIP
        SED_default_filename="/home/rgold/rt/thickdisk7/Te-default/poliresa93.75th130fn6130hi.dat" # WIP
        SED_default=loadtxt(SED_default_filename,usecols=(1,2))

        plot(SED_SCSjet[:,0],SED_SCSjet[:,1],'mx-',label=r"$\rm T_{e,jet}=10m_ec^2,SCS+jet$")
        plot(SED_default[:,0],SED_default[:,1],'cD-',label=r"$\rm default$")

    except IOError:
        pass
        
    # PLOT SED & COMPARE TO OBSERVATIONS
    nu_obs = array([8.45, 14.90, 22.50, 43.00, 87.73, 102., 145., 230.86, 349., 674., 857.]) #, 1500., 3000., 5000.])
    SED_errors = array([0.031, 0.012, 0.015, 0.026, 0.080, 0.1517, 0.2644, 0.1414, 0.1205, 0.3508, 0.2404]) #, 0. , 0., 0.])
    SgrA_SED_Observed=[[],[]] # [nu_Ghz,F_nu] [WIP]
    nu=linspace(0,1000,100)

    def SgrA_SED_FIT(nu):
        return 0.248*nu**0.45 * exp(-(nu/1100.)**2)

    plot(SED_SCS[:,0],SED_SCS[:,1],'bo-',label=r"$\rm T_{e,jet}=35m_ec^2,SCS$")
    errorbar(nu_obs,SgrA_SED_FIT(nu_obs),yerr=SED_errors,fmt='ks',label="observed")
    plot(nu,SgrA_SED_FIT(nu),'k--',alpha=0.5,lw=4,label="FIT") 

    for FILE in FILES_1D:
    
        SED=loadtxt(FILE,usecols=(1,2))
        plot(SED[:,0],SED[:,1],'rd-',label=FILE)# r"$\rm T_{e,jet}=35m_ec^2,SCS$")

    legend(loc="lower right",labelspacing=0.2,fontsize=15)
    xlabel(r"$\nu/Ghz$")
    ylabel(r"$F_\nu/Jy$")
    tight_layout()

##################################
print "============\n====DONE===="
##################################
