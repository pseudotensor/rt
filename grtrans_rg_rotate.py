import grtrans_batch as gr
import pickle
import pylab
from pylab import *
import numpy as np
import copy
import string,scipy,glob,commands
from scipy import ndimage,constants

pc = scipy.constants.parsec # SI
G = scipy.constants.G # SI
c = scipy.constants.c # SI

#t_rg =  # time in secs assuming Sgr A*                                      
Msun = 2e30 # SI
M = 4.3e6 * Msun # SAG A*
rg = G*M/c**2
d_SagA = 8.3e3*pc
rad2microarcsec = 360/(2*pi)*3600*1e6

import astroray_io

## WIP ##
image_size_inM=5. # 10. # in r_s
image_size_rad=image_size_inM * (2*rg)/d_SagA
image_size=image_size_rad * rad2microarcsec

shadow_schwarzschild = sqrt(27)*2*rg # diameter
shadow_maximally_spinning = 9./2.*2*rg # diameter
# shadow_schwarzschild *= image_size
# shadow_maximally_spinning *= image_size
# shadow_schwarzschild *= 0.5 # 20Mx20M
# shadow_maximally_spinning *= 0.5


def plot_shadows(domain):
    '''Plot predicted black hole shadow size in 'xy' or 'uv' domain.'''

    if string.lower(domain)=="xy":
        # gca().add_artist(Circle((npixel/2,npixel/2),radius = shadow_schwarzschild/d_SagA/2. * rad2microarcsec,color="cyan",alpha=0.5,fill=False,lw=4,ls="dashed"))
        gca().add_artist(Circle((npixel/2,npixel/2),radius = shadow_maximally_spinning/2./d_SagA * rad2microarcsec,color="grey",alpha=0.5,fill=False,lw=4,ls="dashed"))
        # gca().add_artist(Circle((0,0),radius = shadow_maximally_spinning/2./d_SagA * rad2microarcsec,color="grey",alpha=0.5,fill=False,lw=4,ls="dashed"))
    elif string.lower(domain)=="uv": # the circles with *radius* uv_schwarzschild in uv correspond to the shadow *diameters*                             
        gca().add_artist(Circle((npixel/2,npixel/2),radius = uv_schwarzschild,color="cyan",alpha=0.5,fill=False,lw=4,ls="dashed"))
        gca().add_artist(Circle((npixel/2,npixel/2),radius = uv_maximally_spinning,color="grey",alpha=0.5,fill=False,lw=4,ls="dashed"))

    return None

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a  
        # simple example...                                                 
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))
# USAGE: imshow(data, norm=MidpointNormalize(midpoint=0), cmap=plt.cm.seismic)                                                                         


xlist=[]
xlist.append(gr.grtrans())
npixel_ASTRORAY=152 # 303 # 152
npixel=100 # 152 # 303 # 152
nr_of_snapshots=1 # 1 # 5 # 100
I_avg=zeros((npixel,npixel))

## SET PARAMETERS FOR GRTRANS ##
nvals=1 # 1: unpolarized, 4: polarization
th=pi/2. # 0. # 0.2 # pi/2. # pi-0.2 # 0.2 # 1.4 # 1.73 # inclination astroray [rad]
# mumin, mumax, nmu -- Minimum/maximum/number of mu = cos(i) of observer camera(s). mu=1,0 corresponds to face-on, edge-on.

# i=0.2rad: 6.5e16~>6.54Jy # ASTRORAY: 6.55Jy
mdot_min=6.5e16 # 1.25e17 # 1e18 # 2e17
mdot_max=mdot_min # 2e17# mdot_min
mdot_nr=1 # 1 # 4
f=857 # 349 # 231
jonfix=1 # 0:off 1:on
extra=1 # 1:on
debug=1 # 1:on

# DEBUGGING: 
# extra=0,debug=1,npixel=1 : Too many iterations ROOT MUST BE BRACKETED FOR ZBRENT 0 10 -1 -Infinity
# extra=0,debug=0,npixel=1 : Too many iterations ROOT MUST BE BRACKETED FOR ZBRENT 0 10 -1 -Infinity
# extra=1,debug=0,npixel=1 : Too many iterations ROOT MUST BE BRACKETED FOR ZBRENT 0 10 -1 -Infinity
# extra=1,debug=1,npixel=1 : Too many iterations ROOT MUST BE BRACKETED FOR ZBRENT 0 10 -1 -Infinity
# extra=1,debug=1,npixel=52: NO OUTPUT, NO ERROR MESSAGE
# extra=0,debug=0,npixel=11 : NO OUTPUT, NO ERROR MESSAGE
# extra=0,debug=0,npixel=50 : WORKS
# extra=1,debug=0,npixel=50 : NO OUTPUT, NO ERROR MESSAGE
# extra=1,debug=1,npixel=50 : NO OUTPUT, NO ERROR MESSAGE
# extra=0,debug=0,npixel=51 : NO OUTPUT, NO ERROR MESSAGE
# extra=0,debug=0,npixel=52 : WORKS

i1=npixel**2/2+int(sys.argv[-1])
i2=i1 # i1+10


I_avg_array = []

# FILE_LIST = glob.glob("astroray-vs-grtrans/shotimag93.75th140f230fn5550case70816925_"+str(npixel_ASTRORAY-1)+"_*view.dat")
# FILE_LIST = glob.glob("astroray-vs-grtrans/shotimag93.75th140f230fn5550case70816925_151_0view.dat")
# FILE_LIST = glob.glob("astroray-vs-grtrans/shotimag93.75th140f230fn5550case70816925_151_4.20973view.dat")
# FILE_LIST = glob.glob("astroray-vs-grtrans/shotimag93.75th140f230fn5550case708169250_151.dat")
# FILE_LIST = glob.glob("astroray-vs-grtrans/shotimag93.75th20f230fn5550case708169251_151.dat")
# FILE_LIST = glob.glob("astroray-vs-grtrans/shotimag93.75th20f857fn5550case708169251_151.dat")
# FILE_LIST = glob.glob("astroray-vs-grtrans/shotimag93.75th20f857fn5550case708169252_151.dat")
# FILE_LIST = glob.glob("astroray-vs-grtrans/shotimag93.75th20f857fn5550case708169253_151.dat") ## <I>~?Jy ## TUNED TO TABULATED Emissivity
# FILE_LIST = glob.glob("astroray-vs-grtrans/shotimag93.75th20f857fn5550case708169254_151.dat") ## <I>~?Jy ## TUNED TO 4th intensity
# FILE_LIST = glob.glob("astroray-vs-grtrans/shotimag93.75th20f857fn5550case708169255_151.dat") ## <I>~?Jy ## TUNED TO 4th intensity
FILE_LIST = glob.glob("astroray-vs-grtrans/shotimag93.75th20f857fn5550case708169256_151.dat") ## <I>~?Jy ## TUNED TO 4th intensity
print "size(FILE_LIST)=",size(FILE_LIST)

IQUV_astroray = zeros((npixel_ASTRORAY,npixel_ASTRORAY,5))
for FILE in FILE_LIST:
    IQUV_astroray += astroray_io.read_astroray_images(FILE,fmin=5550,fmax=5550)/size(FILE_LIST)
    # IQUV_astroray += fliplr(astroray_io.read_astroray_images(FILE,fmin=5550,fmax=5550)/size(FILE_LIST))
    # IQUV_astroray = scipy.ndimage.rotate(IQUV_astroray,-90)

# I_astroray=loadtxt("astroray_image_thickdisk7_4grtrans_comparison.dat")
# I_astroray=loadtxt("astroray_image_thickdisk7_4grtrans_comparison_i20rad.dat").reshape(npixel_ASTRORAY,npixel_ASTRORAY,5)[:,:,0]
# I_astroray=loadtxt(FILE_LIST[0])
# I_astroray=transpose(ndimage.rotate(I_astroray,180-90))
# I_astroray = fliplr(I_astroray)

for mdot in linspace(mdot_min,mdot_max,mdot_nr):

  snapshot=-1
  # I_avg=0. # ?
  I_avg=zeros((npixel,npixel))
  IQUV_avg_grtrans=zeros(([nvals,13+nvals][extra],npixel,npixel))
  # I_avg=empty((npixel,npixel)) # gives apparent agreement

  # for viewing_angle in linspace(0,2.*pi,nr_of_snapshots): # [83:84]:
    # for viewing_angle in arange(0,2.*pi,2.*pi/nr_of_snapshots):
  # for viewing_angle in [4.*pi/3.-pi/2.]: 
  for viewing_angle in [4.*pi/3.]: # not 0.,0.5,0.75
  # for viewing_angle in [-4.*pi/3.]: # not 0.,0.5,0.75
    snapshot+=1

    print "Current viewing_angle =",viewing_angle,"..."

    # For thickdisk data, there is the option to do something Jon suggested a long time ago, which can be turned on with the variable tjonfix=1

    # th + pi/2. ?
    # xlist[-1].write_grtrans_inputs('inputs.in',fname='THICKDISK',nfreq=1,nmu=1,fmin=2.3e11,fmax=2.3e11,ename='POLSYNCHTH',nvals=nvals,spin=0.9375,stansdard=1,nn=[npixel,npixel,400],uout=0.04,mbh=4.3e6, mdotmin=2e+17,mdotmax=2e+17,nmdot=1,mumin=cos(th/180.*pi+pi/2.),mumax=cos(th/180.*pi+pi/2.),gridvals=[-10,10,-10,10],tgfile='dump0000.bin',tdfile='fieldline',tindf=5550,tnt=1,muval=1./101.,toff=0,phi0=viewing_angle)
#     xlist[-1].write_grtrans_inputs('inputs.in',fname='THICKDISK',nfreq=1,nmu=1,fmin=2.3e11,fmax=2.3e11,ename='POLSYNCHTH',nvals=nvals,spin=0.9375,standard=1,nn=[npixel,npixel,400],uout=0.04,mbh=4.3e6, mdotmin=2e+17,mdotmax=2e+17,nmdot=1,mumin=cos(th/180.*pi),mumax=cos(th/180.*pi),gridvals=[-10,10,-10,10],tgfile='dump0000.bin',tdfile='fieldline',tindf=5550,tnt=1,muval=1./101.,toff=0,phi0=viewing_angle)
#    xlist[-1].write_grtrans_inputs('inputs.in',fname='THICKDISK',nfreq=1,nmu=1,fmin=2.3e11,fmax=2.3e11,ename='POLSYNCHTH',nvals=nvals,spin=0.9375,standard=1,nn=[npixel,npixel,400],uout=0.04,mbh=4.3e6, mdotmin=2e+17,mdotmax=2e+17,nmdot=1,mumin=cos(th/180.*pi-pi/2.),mumax=cos(th/180.*pi-pi/2.),gridvals=[-10,10,-10,10],tgfile='dump0000.bin',tdfile='fieldline',tindf=5550,tnt=1,muval=1./101.,toff=0,phi0=viewing_angle)

    # GOOD
    # xlist[-1].write_grtrans_inputs('inputs.in',fname='THICKDISK',nfreq=1,nmu=1,fmin=2.3e11,fmax=2.3e11,ename='POLSYNCHTH',nvals=nvals,spin=0.9375,standard=1,nn=[npixel,npixel,400],uout=0.04,mbh=4.3e6, mdotmin=mdot,mdotmax=mdot,nmdot=1,mumin=cos(th*180./pi),mumax=cos(th*180./pi),gridvals=[-10,10,-10,10],tgfile='dump0000.bin',tdfile='fieldline',tindf=5550,tnt=1,muval=1./101.,toff=0,phi0=viewing_angle)
    # uout->5e-5 <-> rout=20000M
    # BETTER SIZE SCALE
    # xlist[-1].write_grtrans_inputs('inputs.in',fname='THICKDISK',nfreq=1,nmu=1,fmin=2.3e11,fmax=2.3e11,ename='POLSYNCHTH',nvals=nvals,spin=0.9375,standard=1,nn=[npixel,npixel,400],uout=5e-5,mbh=4.3e6, mdotmin=mdot,mdotmax=mdot,nmdot=1,mumin=cos(th*180./pi),mumax=cos(th*180./pi),gridvals=[-10,10,-10,10],tgfile='dump0000.bin',tdfile='fieldline',tindf=5550,tnt=1,muval=1./101.,toff=0,phi0=viewing_angle)
    # xlist[-1].write_grtrans_inputs('inputs.in',fname='THICKDISK',nfreq=1,nmu=1,fmin=2.3e11,fmax=2.3e11,ename='POLSYNCHTH',nvals=nvals,spin=0.9375,standard=1,nn=[npixel,npixel,400],uout=5e-5,mbh=4.3e6, mdotmin=mdot,mdotmax=mdot,nmdot=1,mumin=cos(th*180./pi),mumax=cos(th*180./pi),gridvals=[-10,10,-10,10],tgfile='dump0000.bin',tdfile='fieldline',tindf=5550,tnt=1,muval=1.,toff=0,phi0=viewing_angle,tjonfix=jonfix) # te=tp
    # xlist[-1].write_grtrans_inputs('inputs.in',fname='THICKDISK',nfreq=1,nmu=1,fmin=f*1e9,fmax=f*1e9,ename='POLSYNCHTH',nvals=nvals,spin=0.9375,standard=1,nn=[npixel,npixel,400],uout=5e-5,mbh=4.3e6, mdotmin=mdot,mdotmax=mdot,nmdot=1,mumin=cos(th),mumax=cos(th),gridvals=[-10,10,-10,10],tgfile='dump0000.bin',tdfile='fieldline',tindf=5550,tnt=1,muval=1/101.,toff=0,phi0=viewing_angle,tjonfix=jonfix) # 20160610
    # xlist[-1].write_grtrans_inputs('inputs.in',fname='THICKDISK',nfreq=1,nmu=1,fmin=f*1e9,fmax=f*1e9,ename='POLSYNCHTH',nvals=nvals,spin=0.9375,standard=1,nn=[npixel,npixel,400],uout=5e-5,mbh=4.3e6, mdotmin=mdot,mdotmax=mdot,nmdot=1,mumin=cos(th),mumax=cos(th),gridvals=[-8.6,8.6,-8.6,8.6],tgfile='dump0000.bin',tdfile='fieldline',tindf=5550,tnt=1,muval=1/11.,toff=0,phi0=viewing_angle,tjonfix=jonfix) 
    # xlist[-1].write_grtrans_inputs('inputs.in',fname='THICKDISK',nfreq=1,nmu=1,fmin=f*1e9,fmax=f*1e9,ename='POLSYNCHTH',nvals=nvals,spin=0.9375,standard=1,nn=[npixel,npixel,400],uout=5e-5,mbh=4.3e6, mdotmin=mdot,mdotmax=mdot,nmdot=1,mumin=cos(th),mumax=cos(th),gridvals=[-8.6,8.6,-8.6,8.6],tgfile='dump0000.bin',tdfile='fieldline',tindf=5550,tnt=1,muval=1/11.,toff=0,phi0=viewing_angle,tjonfix=jonfix,debug=1,extra=extra,i1=i1,i2=i2) # debug extra focus on hot pixel 
    # xlist[-1].write_grtrans_inputs('inputs.in',fname='THICKDISK',nfreq=1,nmu=1,fmin=f*1e9,fmax=f*1e9,ename='POLSYNCHTH',nvals=nvals,spin=0.9375,standard=1,nn=[npixel,npixel,400],uout=5e-6,mbh=4.3e6, mdotmin=mdot,mdotmax=mdot,nmdot=1,mumin=cos(th),mumax=cos(th),gridvals=[-8.6,8.6,-8.6,8.6],tgfile='dump0000.bin',tdfile='fieldline',tindf=5550,tnt=1,muval=1/11.,toff=0,phi0=viewing_angle,tjonfix=jonfix,debug=1,extra=extra,i1=i1,i2=i2) # debug extra focus on hot pixel uout 5e-6 ~> check r
    # xlist[-1].write_grtrans_inputs('inputs.in',fname='THICKDISK',nfreq=1,nmu=1,fmin=f*1e9,fmax=f*1e9,ename='POLSYNCHTH',nvals=nvals,spin=0.9375,standard=1,nn=[npixel,npixel,400],uout=1e-3,mbh=4.3e6, mdotmin=mdot,mdotmax=mdot,nmdot=1,mumin=cos(th),mumax=cos(th),gridvals=[-8.6,8.6,-8.6,8.6],tgfile='dump0000.bin',tdfile='fieldline',tindf=5550,tnt=1,muval=1/11.,toff=0,phi0=viewing_angle,tjonfix=jonfix,extra=extra) # uout 1e-3 ~> checked <=r<=668.39419037309119
    # xlist[-1].write_grtrans_inputs('inputs.in',fname='THICKDISK',nfreq=1,nmu=1,fmin=f*1e9,fmax=f*1e9,ename='POLSYNCHTH',nvals=nvals,spin=0.9375,standard=1,nn=[npixel,npixel,400],uout=1e-3,mbh=4.3e6, mdotmin=mdot,mdotmax=mdot,nmdot=1,mumin=cos(th),mumax=cos(th),gridvals=[-8.6,8.6,-8.6,8.6],tgfile='dump0000.bin',tdfile='fieldline',tindf=5550,tnt=1,muval=1/11.,toff=0,phi0=viewing_angle,tjonfix=jonfix,debug=debug,extra=extra,i1=i1,i2=i2) # debug extra focus on hot pixel at 78,52 uout 1e-3 ~> checked <=r<=668.39419037309119

    # xlist[-1].write_grtrans_inputs('inputs.in',fname='THICKDISK',nfreq=1,nmu=1,fmin=f*1e9,fmax=f*1e9,ename='POLSYNCHTH',nvals=nvals,spin=0.9375,standard=1,nn=[npixel,npixel,400],uout=1e-3,mbh=4.3e6, mdotmin=mdot,mdotmax=mdot,nmdot=1,mumin=cos(th),mumax=cos(th),gridvals=[-8.6,8.6,-8.6,8.6],tgfile='dump0000.bin',tdfile='fieldline',tindf=5550,tnt=1,muval=1/11.,toff=0,phi0=viewing_angle,tjonfix=jonfix,debug=debug,extra=extra)
    # xlist[-1].write_grtrans_inputs('inputs.in',fname='THICKDISK',nfreq=1,nmu=1,fmin=f*1e9,fmax=f*1e9,ename='lambda',nvals=nvals,spin=0.9375,standard=1,nn=[npixel,npixel,400],uout=1e-3,mbh=4.3e6, mdotmin=mdot,mdotmax=mdot,nmdot=1,mumin=cos(th),mumax=cos(th),gridvals=[-8.6,8.6,-8.6,8.6],tgfile='dump0000.bin',tdfile='fieldline',tindf=5550,tnt=1,muval=1/11.,toff=0,phi0=viewing_angle,tjonfix=jonfix,debug=debug,extra=extra)
    # xlist[-1].write_grtrans_inputs('inputs.in',fname='THICKDISK',nfreq=1,nmu=1,fmin=f*1e9,fmax=f*1e9,ename='lambda',nvals=nvals,spin=0.9375,standard=1,nn=[npixel,npixel,400],uout=1e-3,mbh=4.3e6, mdotmin=mdot,mdotmax=mdot,nmdot=1,mumin=cos(th),mumax=cos(th),gridvals=[-25,25,-25,25],tgfile='dump0000.bin',tdfile='fieldline',tindf=5550,tnt=1,muval=1/11.,toff=0,phi0=viewing_angle,tjonfix=jonfix,debug=debug,extra=extra)

    xlist[-1].write_grtrans_inputs('inputs.in',fname='HARM',nfreq=1,nmu=1,fmin=2.3e11,fmax=2.3e11,ename='POLSYNCHTH',nvals=1,spin=0.9375,standard=1,nn=[150,150,400],uout=0.04,mbh=4e6, mdotmin=1.57e15,mdotmax=1.57e15,nmdot=1,mumin=.6428,mumax=.6428,gridvals=[-13,13,-13,13],hhfile='dump040',hdfile='dump',hindf=40,hnt=1,muval=1./4.)

    grtrans_file_out = "grtrans_thickdisk7_i"+str(int(th*100))+"rad_view"+str(viewing_angle)+"_mdot"+str(mdot)+"_f"+str(f)+".dat"
    print "grtrans_file_out:",grtrans_file_out

    try:
        IQUV_grtrans=loadtxt(grtrans_file_out).reshape([nvals,13+nvals][extra],npixel,npixel)
    except:
        print "LAUNCHING GRTRANS RUN..."

        xlist[-1].run_grtrans()
        xlist[-1].read_grtrans_output()

        if extra==1:
            IQUV_grtrans=np.transpose(xlist[-1].ivals[:,:,0].reshape((xlist[-1].nx,xlist[-1].ny,nvals+13)))
        else:
            IQUV_grtrans=np.transpose(xlist[-1].ivals[:,:,0].reshape((xlist[-1].nx,xlist[-1].ny,nvals)))

        # commands.getoutput("rm "+grtrans_file_out)
        # savetxt(grtrans_file_out,I)
        savetxt(grtrans_file_out,IQUV_grtrans.flatten())

        # JASON: units of Jy for Sgr A*
        lbh=xlist[-1].inputs.mbh*6.67e-8/9e20*2e33 # MG/c^2*Msun
        F_grtrans = xlist[-1].spec*lbh**2./(8.*3.09e21)**2.*1e23 # 8kpc?

        mdot_code_scale=.0013
        # mdot_cgs = mdot_input * (mdot_code / mdot)
        tcgs=lbh/3e10
        # rhocgs=sp%mdot/mdot/lbh**3*tcgs*rho
        rho=1. # WIP
        # rhocgs=mdot/mdot_code_scale/lbh**3*tcgs*rho

        print "lbh=",lbh,"F_grtrans=",F_grtrans
        # print "rhocgs=",rhocgs
        # lbh= 6.37355555556e+11 F_grtrans= [[ 0.76097975] th=1.4rad mdot=2e17 mu=1./101.
        # ASTRORAY 0.701388540192 Jy (astroray-vs-grtrans/shotimag93.75th140f230fn5550-5550case70816925_151_4.20973view.dat)
        # uout=1e-3 <=> r=668.39419037309119 => F=6.52737856Jy ASTRORAY: F=6.54093164641Jy (0th intensity, 4th intensity is 30% higher
        # rho_norm -> <F>
        # 1e5      -> 5.49Jy (with tabulated emissivity: 4.33Jy)
        # 1.1e5    -> 6.46Jy (with tabulated emissivity: 5.13Jy)
        # 1.105e5  -> 6.51Jy (with tabulated emissivity: 5.17Jy)
        # 1.11e5   -> 6.56Jy (with tabulated emissivity: 5.21Jy)
        # 1.12e5   -> 6.66Jy (with tabulated emissivity: 5.30Jy)
        # 1.2e5    -> 7.49Jy (with tabulated emissivity: 5.98Jy) nu<0 !

    I=IQUV_grtrans[0,:,:]
    


## PLOT IMAGE ## 
    
    BIAS=True
    if BIAS==True:
        plot_idx=1                      # biased
        title_string=["ASTRORAY","GRTRANS","REL DIFF"]
    else:
        plot_idx = random.choice([1,2]) # unbiased
        title_string=["","","REL DIFF"]

    figure(1,figsize=(16.5,4))

    # http://stackoverflow.com/questions/10171618/changing-plot-scale-by-a-factor-in-matplotlib
    # scale=1e-3
    # ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*scale))
    # ax.xaxis.set_major_formatter(ticks)

    clf()
    # subplot(120+plot_idx)
    subplot(130+plot_idx)
    if plot_idx==1:
        plot_idx+=1
    else:
        plot_idx-=1

    print "IQUV_astroray above 99.999th percentile:",IQUV_astroray[IQUV_astroray[:,:,0]>percentile(IQUV_astroray[:,:,0],99.999),0]
    print "I[I>percentile(I,99.999)]",I[I>percentile(I,99.999)]
    outlier_row_idx,outlier_col_idx = where(I>percentile(I,99.999))
    print "OUTLIER AT ROW",outlier_row_idx,"COL",outlier_col_idx

    # CAP=max(amax(I),amax(IQUV_astroray[:,:,0]))
    CAP=max(amax(percentile(I,99.9)),amax(IQUV_astroray[:,:,0]))

    # subplot(131)
    imshow(I,origin='lower',cmap=cm.cubehelix,vmax=CAP) # ,vmax=percentile(I,99))
    cb = colorbar(pad=0)
    cb.formatter.set_scientific(True);cb.formatter.set_powerlimits((0, 0));cb.update_ticks()
    # cb.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    axis('equal')
    axis((0,npixel,0,npixel))
    text(npixel/2,npixel-20,r"$\nu=$"+str(f)+r"$GHz$",fontdict={"color":"yellow","size":15})
    plot_shadows("xy")
    gca().axes.get_xaxis().set_ticks([]);gca().axes.get_yaxis().set_ticks([])
    title(title_string[1])

    # subplot(120+plot_idx)
    subplot(130+plot_idx)

    imshow(IQUV_astroray[:,:,0],origin='lower',cmap=cm.cubehelix,vmax=CAP)
    cb = colorbar(pad=0)
    cb.formatter.set_scientific(True);cb.formatter.set_powerlimits((0, 0));cb.update_ticks()
    axis('equal')
    axis((0,npixel,0,npixel))
    # imshow(I_astroray,origin='lower',cmap=cm.cubehelix)
    text(npixel/2,npixel-20,r"$\nu=$"+str(f)+r"$GHz$",fontdict={"color":"yellow","size":15})
    plot_shadows("xy")
    gca().axes.get_xaxis().set_ticks([]);gca().axes.get_yaxis().set_ticks([])
    title(title_string[0])

    ## approximate emissivities Melrose+
    subplot(133)
    imshow(IQUV_astroray[:,:,4],origin='lower',cmap=cm.cubehelix,vmax=CAP)
    cb = colorbar(pad=0)
    cb.formatter.set_scientific(True);cb.formatter.set_powerlimits((0, 0));cb.update_ticks()
    axis('equal')
    axis((0,npixel,0,npixel))
    # imshow(I_astroray,origin='lower',cmap=cm.cubehelix)
    text(npixel/2,npixel-20,r"$\nu=$"+str(f)+r"$GHz$",fontdict={"color":"yellow","size":15})
    plot_shadows("xy")
    gca().axes.get_xaxis().set_ticks([]);gca().axes.get_yaxis().set_ticks([])
    title(title_string[0]+" (j_approx)")



## ANNOTATIONS ## 
    # tight_layout(pad=0,w_pad=0,h_pad=0)

    I_avg += I/nr_of_snapshots
    IQUV_avg_grtrans += IQUV_grtrans/nr_of_snapshots

    savefig("grtrans_thickdisk7_view"+string.zfill(snapshot,nvals)+"mdot"+str(mdot)+"_f"+str(f)+".png")

  # print "snapshot=",snapshot
  # if (snapshot==nr_of_snapshots-1):
  #   # array of averages for each mdot ...
  #   I_avg_array.append([I_avg])
  #   # I_avg=empty((npixel,npixel)) # gives apparent agreement
  #   # I_avg=zeros((npixel,npixel)) # gives apparent agreement
  #   snapshot=0

I_avg_array = array(I_avg_array[:])



# ###########
# figure(3) #
# ###########

# # title("Averaged images over all viewing angles")
# # clf()
# # imshow(I_avg,origin='lower',cmap=cm.cubehelix)
# X=linspace(-image_size,image_size,npixel)
# Y=linspace(-image_size,image_size,npixel)
# for plot_idx in arange(shape(I_avg_array)[0]):
#     subplot(511+plot_idx)
#     pcolormesh(X,Y,I_avg_array[plot_idx],cmap=cm.cubehelix)
#     axis('scaled')
#     gca().axes.get_xaxis().set_ticks([]);gca().axes.get_yaxis().set_ticks([])
#     plot_shadows("xy")
#     # title(title_string[1])


n_rows=nvals+1 # =5 for polarized case
##########################
# figure(2,figsize=(18,6)) #
# figure(2,figsize=(32/n_rows,10*n_rows)) #
figure(2,figsize=(8,10*n_rows)) # looks ok on GUI but png file looks bad
# figure(2,figsize=(8*3,10*n_rows)) # wastes space in GUI but better as png file :-s
##########################

# title("Averaged images over all viewing angles")
# clf()

# for current_row in range(n_rows):
for current_row in [[0,4],range(n_rows)][nvals>1]:
  # subplot(5,3,1+current_row*3)
  idx=[[0,4],range(n_rows)][nvals>1].index(current_row)
  subplot(n_rows,3,1+idx*3)
  # imshow(I_avg,origin='lower',cmap=cm.cubehelix)

  X=linspace(-image_size,image_size,npixel)
  Y=linspace(-image_size,image_size,npixel)
  # pcolormesh(X,Y,I_avg,cmap=cm.cubehelix)
  CMAP=[cm.cubehelix,cm.PuOr,cm.PuOr,cm.PuOr,cm.cubehelix][current_row]
  #RG:FIXME for nvals lists are inconsistent
  pcolormesh(X,Y,IQUV_avg_grtrans[[0,1,2,3,0][current_row],:,:],cmap=CMAP,vmax=percentile(IQUV_avg_grtrans[[0,1,2,3,0][current_row],:,:],99),norm=[None,MidpointNormalize(midpoint=0),MidpointNormalize(midpoint=0),MidpointNormalize(midpoint=0),None][current_row])
  # colorbar(pad=0,orientation="horizontal")
  plot_shadows("xy")
  plt.axis('off') # gca().axes.get_xaxis().set_ticks([]);gca().axes.get_yaxis().set_ticks([])
  if current_row==0:
      title(title_string[1])
  axis('scaled')

  subplot(n_rows,3,2+idx*3)
  X_astroray = linspace(-image_size,image_size,npixel_ASTRORAY)
  Y_astroray = X_astroray[:]

  # I_avg=amax(I_avg,1e-10)
  # Following is not correct for current_row=4 when I_approximate_astroray should be compared to I_grtrans
  try:
      #RG:FIXME for nvals lists are inconsistent
      # A/<A> - B/<B>
      rel_diff = IQUV_astroray[:,:,current_row]/mean(IQUV_astroray[:,:,current_row])-(IQUV_avg_grtrans[[current_row,0][current_row==4],:,:]/mean(IQUV_avg_grtrans[[current_row,0][current_row==4],:,:])) # /I_avg/mean(I_avg) 
  except:
      pass

  #RG:FIXME for nvals lists are inconsistent
  pcolormesh(X_astroray,Y_astroray,IQUV_astroray[:,:,current_row],cmap=CMAP)
  # colorbar(pad=0,orientation="horizontal")
  plot_shadows("xy")
  plt.axis('off') # gca().axes.get_xaxis().set_ticks([]);gca().axes.get_yaxis().set_ticks([])
  axis('scaled')
  if current_row==0:
      title(title_string[0])

  subplot(n_rows,3,3+idx*3)
  try:
      pcolormesh(X_astroray,Y_astroray,rel_diff,cmap=cm.RdBu,norm=MidpointNormalize(midpoint=0))
  # colorbar(pad=0,orientation="horizontal")
  except:
      pass
  plot_shadows("xy")
  plt.axis('off') # gca().axes.get_xaxis().set_ticks([]);gca().axes.get_yaxis().set_ticks([])
  axis('scaled')
  if current_row==0:
      title(title_string[2])

tight_layout(h_pad=0,w_pad=0,pad=0.0)

# FIXME: DOES NOT SAVE ALL MDOT FILES... :-s
# savefig("grtrans-vs-astroray_thickdisk7_view_avg"+"mdot"+str(mdot)+".png",dpi=400)
savefig("grtrans-vs-astroray_thickdisk7_view_avg"+"mdot"+str(mdot)+".png")
# savefig("grtrans-vs-astroray_thickdisk7_view_avg"+"mdot"+str(mdot)+".pdf")


## IMAGE PROCESSING STUFF TO DETECT SHADOW

# figure(4)
# # http://www.scipy-lectures.org/advanced/image_processing/#edge-detection
# # Use a gradient operator (Sobel) to find high intensity variations:
# # Sobel Filter
# sx = ndimage.sobel(IQUV_astroray[:,:,0], axis=0, mode='constant')
# sy = ndimage.sobel(IQUV_astroray[:,:,1], axis=1, mode='constant')
# sob = np.hypot(sx, sy)
# imshow(sob,cmap=cm.bone)
# title("Sobel filtered (edge detection)")
# axis('off')

if extra==1: # debug==1:

  try:
    import read_geodebug_file as d
  except ImportError:
    pass

  figure(0)
  # d.r d.th d.phi
  x = d.r * sin(d.th) * cos(d.phi)
  y = d.r * sin(d.th) * sin(d.phi)
  z = d.r * cos(d.th)
  tangent_in  = array((diff(x)[ 0],diff(y)[ 0],diff(z)[ 0]))
  tangent_out = array((diff(x)[-1],diff(y)[-1],diff(z)[-1]))
  # A.B=|A||B|cos(angle)
  cos_angle=dot(tangent_in,tangent_out)/norm(tangent_in)/norm(tangent_out)
  print "Deflection angle ",arccos(cos_angle),"rad =",arccos(cos_angle)/2./pi*360.,"deg"
  from mpl_toolkits.mplot3d import Axes3D


  # http://arxiv.org/pdf/gr-qc/9907034v1.pdf eqs.(20),(24)
  r_per=amin(d.r) # units [M]
  b=sqrt(d.alpha**2+d.beta**2) # units [M]
  a=d.a

  deflection_angle_weakfield   = 4./b * ( 1. + 15./16.*pi/b )
  # http://arxiv.org/pdf/1405.2919.pdf eq (3.40) (equatorial), see: Sereno,de Luca for general orbits
  deflection_angle_weakfield_Sereno = 4./b + (15./4./pi-4.*a)/b**2 + (4*a**2-10.*pi*a+128./3.)/b**3 + (15./64.*pi*(76*a**2+231)-4.*a*(a**2+48.))/b**4 + (4.*(a**2+128)*a**2-9./2.*pi*(6.*a**2+77.)*a+3584./5.)/b**5
  deflection_angle_strongfield = log( 3.482/(b-3.*sqrt(3.)))
  # https://arxiv.org/pdf/gr-qc/0611086v2.pdf eqs.(8),(14->23),(17->25)

  alpha_vs_b_file = file("grtrans_alpha_vs_b.dat","a")
  savetxt(alpha_vs_b_file,array([i1,i2,d.alpha,d.beta,b,arccos(cos_angle),deflection_angle_weakfield_Sereno]).T,newline=" ")
  # savetxt(alpha_vs_b_file,["\n"])
  alpha_vs_b_file.write("\n")
  alpha_vs_b_file.close()

  # http://arxiv.org/pdf/1405.2919.pdf see Figs. 2.2,3.5
  # ASTRORAY_geodesic.dat
  print "pericenter=",r_per
  print "impact parameter b=",b
  print "...assuming Schwarzschild..."
  print "deflection_angle_weakfield=",deflection_angle_weakfield,"rad =",deflection_angle_weakfield/2./pi*360.,"deg"
  print "deflection_angle_weakfield_Sereno=",deflection_angle_weakfield_Sereno,"rad =",deflection_angle_weakfield_Sereno/2./pi*360.,"deg"
  print "deflection_angle_strongfield=",deflection_angle_strongfield,"rad =",deflection_angle_strongfield/2./pi*360.,"deg"

  figure(-1,figsize=(8,14))
  subplot(511)
  ticklabel_format(style="sci",scilimits=(14,1),useOffset=True)
  plot(d.lam,d.ii,'.') # affine parameter vs intensity ?
  xlabel("affine parameter (?)") # affine parameter vs intensity ?
  ylabel("total intensity (?)") # affine parameter vs intensity ?

  subplot(512)
  ticklabel_format(style="sci",scilimits=(14,1),useOffset=True)
  plot(d.lam,d.ji,'.') # affine parameter vs j_nu_I ?
  xlabel("affine parameter (?)") # affine parameter vs intensity ?
  ylabel(r"$j_{\nu,I}$") 

  subplot(513)
  ticklabel_format(style="sci",scilimits=(14,1),useOffset=True)
  plot(d.r,d.ii,'.') # radius vs intensity ?
  xlabel(r"$r/M$") 
  ylabel("total intensity (?)") # intensity ?

  subplot(514)
  ticklabel_format(style="sci",scilimits=(14,1),useOffset=True)
  plot(d.r,d.taui,'.') # radius vs optical depth ?
  xlabel(r"$r/M$") 
  ylabel(r"$\tau_I$")

  subplot(515)
  ticklabel_format(style="sci",scilimits=(14,1),useOffset=True)
  plot(d.r,d.rho,'.-') # radius vs rest-mass density ?
  xlabel(r"$r/M$") 
  ylabel(r"$\rho$")


  tight_layout(pad=0)

  print "d.mdot=",d.mdot


