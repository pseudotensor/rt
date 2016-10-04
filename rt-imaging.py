# VISUALIZE IMAGES OBTAINED BY THE ASTRORAY CODE
# Roman Gold 2014

#################################################
# Broderick et al http://arxiv.org/abs/1311.5564:
# 
# Sgr A*: 
# mass (4.3 /pm 0.5) x 10^6 Msun 
# distance (8.3 \pm 0.4) kpc
# lensed image of Sgr A* horizon is 53 \pm 2 microarcs (Ghez et al. 2008; Gillessen et al. 2009a,b)
# resolution of existing and future mm-VLBI: (30--10) microarcsec 
#
# Bardeen 1973, Luminet 1979, Johannsen & Psaltis 2010
# Shadow diameter:
# a=0: sqrt(27) R_s
# a=1: 9/2 R_s
#########################################################

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

rc('font',size=18)

################
## USER SPECS ##
# TIME_AVERAGE=True # generalized code, obsolete now
OUTPUT_EXT="png" # "pdf" # A WORD OF CAUTION: "pdf" option is very slow...
miniversion = True # False True
WANTED_PLOTS=["IQUV","IP"] # "IP","IQUV"
# WANTED_PLOTS=["IP"]

filename = sys.argv[1] # should point to shotimage*.dat file

HOME=commands.getoutput("echo $HOME")
# RT_DIR="/rt/"
# RT_DIR="/codes/rt-git/"
RT_DIR="/codes/astroray/"

I_xy_scale=1
Jy2cgs=1e23
ang_size_norm=66.4648/Jy2cgs

# assumes obs.txt (or similar file) in same dir (provided by Andrew Chael see [eht_python_for_roman.zip])
# SOMEWHAT HARDWIRED DIRECTORY...
try:
    EHT_config_file="EHT-obs-"+filename.split('fn')[0].split('f')[1]+"Ghz-2017.txt"
    # EHT_config_file="EHT-obs-230Ghz-2017.txt"
    # EHT_config_file="obs.txt"
    # EHT_config_file="obs-SMT-SMA.txt"
    eht_obs_uv = loadtxt(HOME+RT_DIR+EHT_config_file,usecols=[0,4,5],comments='#')
    print "SUCCESSFULLY READ IN EHT uv-TRACKS"
except:
    try:
        RT_DIR="/rt/"
        eht_obs_uv = loadtxt(HOME+RT_DIR+EHT_config_file,usecols=[0,4,5],comments='#')
        print "SUCCESSFULLY READ IN EHT uv-TRACKS"
    except:
        print "Could not load EHT array configuration file! Ignore..."
        pass



## Taken from Michael's Notebook ##
SGRAdec = -(29 + 28.118/3600); ## (* Source Parameters *)
Lambda = 1.3e-3;            ## (* Wavelength (meters) *)

# (* Telescope Parameters *)
SMAvec = array([-5464523.400, -2493147.080, 2150611.750])/Lambda
SMTvec = array([-1828796.200, -5054406.800, 3427865.200])/Lambda
CARMAvec = array([-2397431.300, -4482018.900, 3843524.500])/Lambda
LMTvec = array([-768713.9637, -5988541.7982, 2063275.9472])/Lambda
ALMAvec = array([2225037.1851, -5441199.1620, -2479303.4629])/Lambda
PVvec = array([5088967.9000, -301681.6000, 3825015.8000])/Lambda
PdBIvec = array([4523998.40 , 468045.240 , 4460309.760])/Lambda
SPTvec = array([0.0 , 0.0, -6359587.3])/Lambda
Haystackvec = array([1492460, -4457299, 4296835])/Lambda
#######################################################
# SGRAvec = [cos[SGRAdec*pi/180], 0, sin[SGRAdec*pi/180]]
# projU = Cross[[0, 0, 1], SGRAvec];
# projU = projU/Norm[projU];
# projV = -Cross[projU, SGRAvec];
# RR[vec_, \[Theta]_] := [vec[[1]]*Cos[\[Theta]] - 
#     vec[[2]]*Sin[\[Theta]], 
#    vec[[1]]*Sin[\[Theta]] + vec[[2]]*Cos[\[Theta]], vec[[3]]];
# RRelevcut[vec_, \[Theta]_] := 
#  If[RR[vec, \[Theta]].SGRAvec/Norm[vec]/Norm[SGRAvec] > 
#    Cos[75*Pi/180], 1, 
#   0.0]; (* 15 Degree Elevation Cut *)
####################################################
# Cross[projU, projV] ## (* Ensure that the coordinate system is right handed *)
# SGRAvec
# {0.8745536098088065, 0, -0.4849288438218393}
# {0.874554, 0, -0.484929}
# {0.874554, 0., -0.484929}
####################################################


if "vary" in sys.argv: # use that for movies where parameter varies (adds current parameter value to plot titles)
    VARY = string.lower(sys.argv[sys.argv.index("vary")+1])
else:
    VARY="" # can be "","magn_cap","theta","r"

SCATTERING = "ON"
try: 
    SCATTERING_REDUCTION=float(sys.argv[sys.argv.index("blur")+1])
    print "SCALING SCATTERING KERNEL (BLURRING) BY ",SCATTERING_REDUCTION
except:
    SCATTERING_REDUCTION=0.5 # 1: full scattering
if SCATTERING == "ON":
    try:
        from scipy import ndimage
    except:
        pass

# help(ndimage.filters)
# FUNCTIONS
#     convolve(input, weights, output=None, mode='reflect', cval=0.0, origin=0)
#         Multi-dimensional convolution.
        
#         The array is convolved with the given kernel.

# see for Gaussian and other Kernels:
# http://nbviewer.ipython.org/github/agile-geoscience/notebooks/blob/master/Filtering_horizons.ipynb
# http://nbviewer.ipython.org/github/mroberts3000/GpuComputing/blob/master/IPython/GaussianBlur.ipynb

# DETECT EDGE OF SHADOW SIMILAR TO DIMITRIOS USING:
# ndimage.gaussian_gradient_magnitude

# great image detection tutorial:
# http://pythonvision.org/basic-tutorial/
# http://scipy-lectures.github.io/advanced/image_processing/

angle_unit="arcsec" # "rad"

# This is used for enhancing smoother uv data at large scales. 
# BE CAREFUL when the xy-data don't decay to zero for large xy this introduces artefacts
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
shadow_schwarzschild = sqrt(27)*2*rg # diameter
shadow_maximally_spinning = 9./2.*2*rg # diameter
#######################################################


## SETTINGS ##
#nxy = 201 # image resolution (nxy+1)x(nxy+1) see ASTRO_RAY_main.cpp
iter = filename.split("fn")[1].split("case")[0]

## READ-IN ##
fp = open(filename,"rb")
header = fromfile(fp,count=20)
fp.close()
nxy=int(header[2])
observing_frequency=header[3]
image_size_inM=header[4]
image_size_rad=image_size_inM * (2*rg)/d_SagA
image_size=image_size_rad * rad2microarcsec

pixeldim = image_size/nxy # Specify the linear size of a pixel, in \[Mu]as
if angle_unit=="rad":
    pixeldim = image_size_rad/nxy

# WIP: UNDERSTAND THIS!
# bh:
X = pixeldim*arange(-round(nxy/2)-1,round(nxy/2)+1)
# laptop:
# X = pixeldim*arange(-round(nxy/2),round(nxy/2)+1)

#?X = pixeldim*arange(-round(nxy/2),round(nxy/2))
Y = X[:]

freq_unit=1e-9 # uv plane scale
uvspacing = image_size_rad/nxy


data = empty((nxy+1,nxy+1,5))
IMAGE_FILES = [entry for entry in sys.argv[1:] if "shotimag" in entry]
filename=filename.replace("fn"+iter,"fn"+iter+"-"+IMAGE_FILES[-1].split("fn")[1].split("case")[0])
print filename
for filename_snapshot in IMAGE_FILES:
    fp = open(filename_snapshot,"rb")
    header = fromfile(fp,count=20)
    nxy=int(header[2])
    data += fromfile(fp,dtype=float64).reshape(nxy+1,nxy+1,5)/size(IMAGE_FILES)
    fp.close()
filename_out = filename.replace(".dat","."+OUTPUT_EXT)
# filename_out = filename.replace(".dat",".png")

## HEADER INFO (see [imaging.cpp]) ##
a=header[0];th=header[1];nxy=int(header[2]);
 #    header[3]=double(sftab[kk][0]);
 #    header[4]=double(sftab[kk][1]);
heat=header[5];rhonor=header[6]
I_img_avg=header[7];LP_img_avg=header[8];EVPA_img_avg=header[9];CP_img_avg=header[10]
 #    header[11]=double(err[kk]);
TpTe=header[12]
mdot=header[13] # in [year/Msun]

print "Image-averaged (zero-baseline) flux: ",I_img_avg,"Jy"
print "Image (pixel) resolution: ",pixeldim,"muas"

try:
    if VARY=="r":
        rr = 1+( float(filename_out.split(str(int(nxy))+"_")[1].split('.')[0]) -100)/5.
        title_vary_string = " at $r="+str(rr)+"M$"
    elif VARY=="t":
        iter_ref=2010 # 6100
        dt=5 # 4
        t = ( float(filename_out.split("_")[0].split('fn')[1]) - iter_ref )*dt *rg/c /60.
        title_vary_string = " at $t="+str(int(round(t,0)))+"min$"
    elif VARY=="theta":
        th = pi*(float(filename_out.split(str(int(nxy))+"_")[1].split('.')[0])-100)/100.
        title_vary_string = " at $\Theta="+str(around(th/pi,1))+"\pi$"
    elif VARY=="magn_cap":
        magn_cap = 1.+(float(filename_out.split(str(int(nxy))+"_")[1].split('.')[0])-100)/10.
        title_vary_string = r" at $b^2/\rho <"+str(around(magn_cap,1))+"$"
    elif VARY=="magn_floor":
        magn_floor = (float(filename_out.split(str(int(nxy))+"_")[1].split('.')[0])-100)/100.
        title_vary_string = r" at $b^2/\rho >"+str(around(magn_floor,2))+"$"
    else:
        print "DID NOT UNDERSTAND WHAT PARAMETER YOU ARE VARYING. DISABLE INFO."
        VARY,title_vary_string="",""

except:
    title_vary_string = ""

try: 
    # TRY cm.cubehelix (incremental)
    # Green, D. A., 2011, `A colour scheme for the display of 
    # astronomical intensity images', 
    # Bulletin of the Astronomical Society of India, 39, 289.
    # (2011BASI...39..289G at ADS.) 
    # Please cite this paper if you use `cubehelix' in any publications.
    ####################################################################
    # For diverging colormaps see:
    # http://dx.doi.org/10.1007/978-3-642-10520-3_9
    # and notice issue of shifted, diverging colormap
    # http://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib

    #colormaps = [cm.gnuplot2,cm.PuOr,cm.bwr,cm.RdBu_r]
    # QUANTITY  [I           ,Q      ,U      ,V      ]
    colormaps = [cm.cubehelix,cm.PuOr,cm.PuOr,cm.PuOr]
    colormaps_IQUV_uv = [cm.cubehelix,cm.gnuplot2,cm.gnuplot2,cm.gnuplot2]
    #colormaps_4panel = [cm.cubehelix,cm.gnuplot2,cm.RdBu_r,cm.bone]
    colormaps_miniversion_4panel = [cm.cubehelix,cm.afmhot,cm.RdBu_r,cm.bone]
    colormaps_4panel = [cm.cubehelix,cm.gnuplot2,cm.RdBu_r,cm.PuOr]
    #colormaps_miniversion_4panel = [cm.cubehelix,cm.afmhot,cm.RdBu_r,cm.PuOr]
except AttributeError: # only want to catch unavailable colormaps
    colormaps = [cm.afmhot,cm.afmhot,cm.afmhot,cm.afmhot]
    colormaps_4panel = [cm.afmhot,cm.afmhot,cm.afmhot,cm.afmhot]
    colormaps_miniversion_4panel = [cm.afmhot,cm.afmhot,cm.afmhot,cm.afmhot]

# CREDIT to Joe Kington's answer in:
# http://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib
 
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

def fmt(x, pos):
    '''Customized formatter. Credit goes to the answer in:
    http://stackoverflow.com/questions/25983218/scientific-notation-colorbar-in-matplotlib'''
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${}\times 10^{{{}}}$'.format(a,b) # how to save the white space between a and \times ?
    #return r'${}\times 10^{{{}}}$'.format(a,b)

# control figure placement on screen 
manager = get_current_fig_manager()
fig_pos=["+0+0","+500+0","+0+500","+500+500","+250+250"] # Need to understand syntax better... 
titles = ["I","Q","U","V"]
titles_4panel_xy = [ 
r"$\rm I$", # \times 10^"+str(int(log10(I_xy_scale)))+"$", 
r"$\sqrt{|Q|^2+|U|^2}$", # \times 10^"+str(int(log10(I_xy_scale)))+"$", 
r"$\rm EVPA$", # =\pi/2. - 0.5*angle(Q_{xy}+iU_{xy})$", ## TESTED
r"$\rm V$"] # \times 10^"+str(int(log10(I_xy_scale)))+"$"
titles_IQUV_uv = [ 
r"$\rm |\tilde{I}|$", 
r"$\rm |\tilde{Q}|$",
r"$\rm |\tilde{U}|$",
r"$\rm |\tilde{V}|$"]
titles_4panel_uv = [ 
# r"$\rm \breve{PF}\equiv |\tilde{Q}+i\tilde{U}|$", 
r"$\rm |\tilde{I}|$", 
r"$\rm |\breve{m}|\equiv |(\tilde{Q}+i\tilde{U}) / \tilde{I}|$", 
r"$\rm EVPA$", # \equiv phase(\tilde{P}/\tilde{I}) \times 90/\pi$", 
r"$\rm |\breve{v}|\equiv |\tilde{V} / \tilde{I}|$"]
titles_miniversion_4panel_uv = [ 
r"$|\rm \tilde{I}|$", 
r"$\rm |\breve{m}|\equiv | (\tilde{Q}+i\tilde{U}) / \tilde{I}|$", 
r"$\rm EVPA$", # \equiv phase(\tilde{P}/\tilde{I}) \times 90/\pi$", 
r"$\rm |\breve{v}|\equiv |\tilde{V} / \tilde{I}|$"]
# r"$|\tilde{\rm m}_{LP}|\equiv|\mathcal{FFT}\{(Q+iU)/I\}|$"]

# 150microarcsec = 1.4 in u-v plane 
# 15 microarcsec = 14  in u-v plane
# for plot labels given the
# nominal scaling of 10^9\lambda for the baseline: shadow size .ie. diameter
uv_schwarzschild = 150.*1.4/(shadow_schwarzschild /d_SagA * rad2microarcsec)
uv_maximally_spinning = 150.*1.4/(shadow_maximally_spinning /d_SagA * rad2microarcsec)



def plot_shadows(domain):
    '''Plot predicted black hole shadow size in 'xy' or 'uv' domain.'''

    if string.lower(domain)=="xy":
        gca().add_artist(Circle((0,0),radius = shadow_schwarzschild/d_SagA/2. * rad2microarcsec,color="cyan",alpha=0.5,fill=False,lw=4,ls="dashed"))
        gca().add_artist(Circle((0,0),radius = shadow_maximally_spinning/2./d_SagA * rad2microarcsec,color="grey",alpha=0.5,fill=False,lw=4,ls="dashed"))
    elif string.lower(domain)=="uv": # the circles with *radius* uv_schwarzschild in uv correspond to the shadow *diameters*
        gca().add_artist(Circle((0,0),radius = uv_schwarzschild,color="cyan",alpha=0.5,fill=False,lw=4,ls="dashed"))
        gca().add_artist(Circle((0,0),radius = uv_maximally_spinning,color="grey",alpha=0.5,fill=False,lw=4,ls="dashed"))

    return None

def plot_EHT_uv_tracks(config='2017'):
    scatter(eht_obs_uv[:,1]/1e9,eht_obs_uv[:,2]/1e9,s=10,c="r",marker="o",alpha=0.25,label="EHT 2017")
    return None


def plot_polticks(every=6,scale_quiver=1e-5,width_quiver=0.01,I_threshold=0.05):
    '''PLOT polarization ticks indicating magnitude and direction of linear polarization'''

    ## WIP ##
    # I mask?, edges around ticks (arrow shafts)?
    # every=5, I_threshold=0.1, scale_quiver=0.2, width_quiver=0.002

    ## WORKS WELL FOR 230Ghz (SGR A*)
    # every=6
    # I_threshold=0.05 
    # scale_quiver=1e-4 # 0.2 * amax(I_xy)  # 0.2  good for LP/I
    # width_quiver=0.01 # 0.01 good for LP/I
    ## WORKS WELL FOR 102Ghz (SGR A*)
    # every=6
    # I_threshold=0.05 
    # scale_quiver=1e-2*amax(I_xy)  # 0.2  good for LP/I
    # width_quiver=0.01 # 0.01 good for LP/I

    # scale_quiver=1e-5   # M87 & AVERY's MODEL

    # quiver_fake_inst= quiver(X[::every],Y[::every],Q_masked[::every,::every],U_masked[::every,::every],headlength=0.,headaxislength=0.,headwidth=0.,color="silver",alpha=0.9,pivot="mid",width=width_quiver/scale_quiver,scale_units="xy",angles="uv",scale=scale_quiver) # ,angles=EVPA_xy[::every,::every])
    # NIRVANA=1000

    # DEFAULT
    # quiver_fake_inst= quiver(X[::every]+NIRVANA,Y[::every]+NIRVANA,(Q_xy/I_xy)[::every,::every],(U_xy/I_xy)[::every,::every],headlength=0.,headaxislength=0.,headwidth=0.,color="silver",alpha=0.9,pivot="mid",width=width_quiver/scale_quiver,scale_units="xy",angles="xy",scale=scale_quiver) # ,angles=EVPA_xy[::every,::every])

    # quiver_fake_inst= quiver(X[::every]+NIRVANA,Y[::every]+NIRVANA,(Q_xy/I_xy)[::every,::every],(U_xy/I_xy)[::every,::every],headlength=0.,headaxislength=0.,headwidth=0.,color="silver",alpha=0.9,pivot="mid",scale_units="xy",angles="xy") # ,angles=EVPA_xy[::every,::every])

    EVPA_xy = pi/2. - 0.5*angle(Q_xy+1j*U_xy)
    x_michael = abs(Q_xy+1j*U_xy)*cos( pi/2 - 0.5*angle(Q_xy+1j*U_xy) )
    y_michael = abs(Q_xy+1j*U_xy)*sin( pi/2 - 0.5*angle(Q_xy+1j*U_xy) )

    # x_michael/=I_xy;y_michael/=I_xy

    Q_masked=x_michael
    Q_masked[I_xy < I_threshold*amax(I_xy)] = None # 0.
    U_masked=y_michael
    U_masked[I_xy < I_threshold*amax(I_xy)] = None # 0.

    # silver color
    # quiver_inst2 = quiver(X[::every],Y[::every],Q_masked[::every,::every],U_masked[::every,::every],headlength=0.,headaxislength=0.,headwidth=0.,color="silver",alpha=0.5,pivot="mid",width=width_quiver,scale=scale_quiver,angles="xy")
    # quiver_inst2 = quiver(X[::every],Y[::every],Q_masked[::every,::every],U_masked[::every,::every],headlength=0.,headaxislength=0.,headwidth=0.,color="black",alpha=0.5,pivot="mid",width=width_quiver,scale=scale_quiver,angles="xy")
    quiver_inst2 = quiver(X[::every],Y[::every],Q_masked[::every,::every],U_masked[::every,::every],headlength=0.,headaxislength=0.,headwidth=0.,color="y",alpha=0.5,pivot="mid",width=width_quiver,scale=scale_quiver,angles="xy")
    # quiver_inst2 = quiver(X[::every],Y[::every],Q_masked[::every,::every],U_masked[::every,::every],headlength=0.,headaxislength=0.,headwidth=0.,color="orange",alpha=0.5,pivot="mid",width=width_quiver,scale=scale_quiver,angles="xy")

    # quiverkey(quiver_inst2, 0.7, 0.92, 0.5, r'$50\%$', coordinates='figure', labelpos='W',fontproperties={'weight': 'bold', 'size': 15})

    # quiverkey(quiver_inst2, 0.85, 0.97, amax(LP_xy), r'$LP='+str(round(amax(LP_xy)*LP_xy_scale,2))+'Jy/\mu as$', coordinates='figure', labelpos='W',fontproperties={'weight': 'bold', 'size': 15},labelcolor='grey')
    # above axes (upper left)
    # quiverkey(quiver_inst2, 0.25, 1.05, amax(LP_xy), r'$LP:'+str(round(amax(LP_xy)*LP_xy_scale,2))+r'Jy/\mu as^2}$', coordinates='axes', labelpos='W',fontproperties={'weight': 'bold', 'size': 12},labelcolor='y')
    quiverkey(quiver_inst2, -42, 36, amax(LP_xy), r'$LP:'+str(round(amax(LP_xy)*LP_xy_scale,2))+r'Jy/\mu as^2}$', coordinates='data', labelpos='E',fontproperties={'weight': 'bold', 'size': 15},labelcolor='y')

    # quiverkey(quiver_inst, 0.8, 0.98, 0.5, r'$50\%$', coordinates='figure', labelpos='W',fontproperties={'weight': 'bold', 'size': 15})


    return None


##############
# PROCESSING #
##############

if string.lower(SCATTERING)=="on" or string.lower(SCATTERING)=="yes":
    # r_fresnel = # 300 km [Goodman & Narayan 1989]
    # fresnel_kernel = exp(1j*(r/(2.*r_fresnel))**2.)

    # appropriate value for 1.3mm SagA*
    # See these references in Fish et al 1409.4690:
    # Bower et al (2006): 1.039 x 0.64 marcsec
    # Shen  et al (2005): 1.41  x 0.75 marcsec
    # @1.3mm point source is blurred to 22muarcsec
    # http://www.cfa.harvard.edu/sma/events/smaConf/posters/images/Doeleman_SMA10_talk.pdf (slide 11 on piercing scattering screen).
    # http://adsabs.harvard.edu/abs/2009A&A...496...77F (Figure 1)
    # http://astro.berkeley.edu/~gbower/ps/sgravlba2.pdf (section on scattering)
    # FWHM(lambda=1cm) = 1mas and FWHM = 2.3*sigma
    # sigma_in_microarcsec = 7.4 # valid for f=230Ghz
    FWHM = 1e3*((c/observing_frequency*freq_unit)/1e-2)**2. # valid for general frequency
    sigma_in_microarcsec = FWHM/2.3
    sigma_in_microarcsec *= SCATTERING_REDUCTION # Scattering can be partially undone, so... REF: http://adsabs.harvard.edu/abs/2014ApJ...795..134F
    sigma = sigma_in_microarcsec / (image_size/nxy) # 10. # pixel units for ndimage.gaussian_filter

    #XY_2D = meshgrid(X,Y) # shape=(2, nxy, nxy)
    # "kernel = exp(-(X**2+Y**2)/2./sigma)" ## Gaussian
    #kernel = exp(-(XY_2D[0]**2+XY_2D[1]**2)/2./sigma) ## Gaussian

###########################
# http://nbviewer.ipython.org/github/agile-geoscience/notebooks/blob/master/Filtering_horizons.ipynb
    # import scipy.signal
    # new_output = scipy.signal.convolve2d(noisy_horizon, kernel)
    # plt.imshow(new_output, aspect=0.5, vmin=vmin, vmax=vmax)
############################

    I_xy = data[:,:,0]
    for CHANNEL in range(shape(data)[-1]):
        data[:,:,CHANNEL] = ndimage.gaussian_filter(data[:,:,CHANNEL],sigma=sigma) ## sigma is in pixel units
    #    CHANNEL *= kernel
    #    CHANNEL *= fresnel_kernel
        # 
        # data[:,:,CHANNEL]=data[:,:,CHANNEL]/mean(I_xy)*I_img_avg # get I,Q,U,V units in Jy, such that mean(I_xy)=I_img_avg

units_Jy=I_img_avg/mean(data[:,:,0])

for i in range(4):
    data[:,:,i] *= units_Jy/pixeldim**2
    data[:,:,i] /= nxy*nxy

I_xy = data[:,:,0]
Q_xy = data[:,:,1]
U_xy = data[:,:,2]
V_xy = data[:,:,3]

# EVPA_xy = arctan2(Q_xy,U_xy)*90./pi
EVPA_xy = 180./2. - 0.5*angle(Q_xy+1j*U_xy,deg=True)

mbreve_xy = abs((Q_xy+1j*U_xy)/I_xy)
m_xy = (Q_xy+1j*U_xy)/I_xy
# I_xy_masked = I_xy[:,:]
# I_xy_masked[I_threshold_mask]=0

#I_threshold_mask = I_xy >= 0. # NO MASK...
#I_threshold_mask = I_xy < 2e-4 # 0.1*average(I_xy)

# print "Artifically lowering polarization. See whether mbreve is still high."
# Q_xy /= 5.
# U_xy /= 5.
# V_xy /= 5.
# Q_xy[I_threshold_mask] = 0.
# U_xy[I_threshold_mask] = 0.
# V_xy[I_threshold_mask] = 0.

# Q_xy = clip(Q_xy,0.1*average(I_xy),None)
# U_xy = clip(Q_xy,0.1*average(I_xy),None)
# V_xy = clip(Q_xy,0.1*average(I_xy),None)

def w(arr):
    '''Apply blackman window function along both axis of 2d array (removes artefacts due to non-periodic boundaries).'''
    try:
        arr = mlab.apply_window(arr,blackman(shape(arr)[0]),axis=0)
        out = mlab.apply_window(arr,blackman(shape(arr)[1]),axis=1)
    except:
        out = arr
    return out

# I_xy_mean=mean(I_xy);Q_xy_mean=mean(Q_xy);U_xy_mean=mean(U_xy);V_xy_mean=mean(V_xy)
I_uv = fftpack.fftshift(fftpack.fft2(w(I_xy),shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
Q_uv = fftpack.fftshift(fftpack.fft2(w(Q_xy),shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
U_uv = fftpack.fftshift(fftpack.fft2(w(U_xy),shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
V_uv = fftpack.fftshift(fftpack.fft2(w(V_xy),shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))

# ORIGINAL:
# I_uv = fftpack.fftshift(fftpack.fft2(I_xy,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
# Q_uv = fftpack.fftshift(fftpack.fft2(Q_xy,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
# U_uv = fftpack.fftshift(fftpack.fft2(U_xy,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
# V_uv = fftpack.fftshift(fftpack.fft2(V_xy,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))

# u = fftfreq(shape(I_uv)[0],d=uvspacing)*freq_unit
# v = fftfreq(shape(I_uv)[1],d=uvspacing)*freq_unit
# u_incr = unique(fftfreq(shape(I_uv)[0],d=uvspacing)*freq_unit)
# v_incr = unique(fftfreq(shape(I_uv)[1],d=uvspacing)*freq_unit)
u = unique(fftfreq(shape(I_uv)[0],d=uvspacing)*freq_unit)
v = unique(fftfreq(shape(I_uv)[1],d=uvspacing)*freq_unit)


mbreve_uv = abs((Q_uv+1j*U_uv)/I_uv)
P_uv = fftpack.fftshift(fftpack.fft2(sqrt(Q_xy**2+U_xy**2),shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
LP_xy = abs(Q_xy+1j*U_xy)
Pbreve_uv = Q_uv+1j*U_uv
mtilde_uv = fftpack.fftshift(fftpack.fft2((Q_xy+1j*U_xy)/I_xy,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
mtilde_LP_uv = fftpack.fftshift(fftpack.fft2((Q_xy+1j*U_xy)/I_xy))

u_no_zeropadding = unique(fftfreq(shape(mtilde_LP_uv)[0],d=uvspacing)*freq_unit)
v_no_zeropadding = unique(fftfreq(shape(mtilde_LP_uv)[1],d=uvspacing)*freq_unit)

# find the mag and phase
I_mag = abs(I_uv)
I_phase = angle(I_uv)
## WRONG: EVPA_uv = angle(P_uv/I_uv) * 90./pi ## Michael: EVPA
EVPA_uv = angle(Pbreve_uv/I_uv) * 90./pi ## Michael: EVPA CORRECT

# GLOBAL intensity scale to save space with tick labels
# I_xy_scale=1 # 1e4
I_xy_scale=10**-round(log10(amax(I_xy))) # automatic nearest power 10
I_uv_scale=10**-round(log10(amax(abs(I_uv)))) # automatic nearest power 10
Q_xy_scale=10**-round(log10(amax(Q_xy))) # automatic nearest power 10
U_xy_scale=10**-round(log10(amax(U_xy))) # automatic nearest power 10
LP_xy_scale=10**-round(log10(amax(LP_xy))) # automatic nearest power 10
V_xy_scale=10**abs(round(log10(amax(V_xy)))) # automatic nearest power 10
# Not good... amin(V_xy)<0 !
# V_xy_scale=max(V_xy_scale,10**abs(round(log10(amin(V_xy))))) # automatic nearest power 10
V_xy_scale=max(V_xy_scale,10**round(log10(amax(abs(V_xy))))) # automatic nearest power 10


## PRODUCE ARTIFICIAL RING DATA IN uv SPACE ##
ring_data_fourier_space = empty((size(u_no_zeropadding),size(v_no_zeropadding)))
for i in range(shape(u_no_zeropadding)[0]):
    for j in range(shape(v_no_zeropadding)[0]):
        R = sqrt(u_no_zeropadding[i]**2+v_no_zeropadding[j]**2)
        if R>=3. and R<=5.:
            ring_data_fourier_space[i][j] += 1
        else:
            ring_data_fourier_space[i][j] = 0
##############################################OA


###########
## PLOTS ##
###########


# RG: None does not work because ticks are set based on these values
# RG: FINISH data depend limits
# figure(0):             [I_xy          ,Q_xy       ,U_xy       ,V_xy      ]
limits_colors_xy = array([(0,amax(I_xy)*I_xy_scale),(amin(Q_xy)*Q_xy_scale,amax(Q_xy)*Q_xy_scale),(amin(U_xy)*U_xy_scale,amax(U_xy)*U_xy_scale),(amin(V_xy)*V_xy_scale,amax(V_xy)*V_xy_scale)])
# limits_colors_xy = array([(0,amax(I_xy)*I_xy_scale),(amin(Q_xy)*Q_xy_scale,amax(Q_xy)*Q_xy_scale),(amin(U_xy)*U_xy_scale,amax(U_xy)*U_xy_scale),(amin(V_xy)*V_xy_scale,amax(V_xy)*V_xy_scale)])
# limits_colors_xy = [(0,8e-4),(amin(Q_xy),amax(Q_xy)),(amin(U_xy),amax(U_xy)),(amin(V_xy),amax(V_xy))]
#limits_colors_xy = [(0,6e-4),(-4e-5,4e-5),(-4e-5,4e-5),(-1e-5,1e-5)]
# # figure(8):              [I_xy    ,mbreve_xy,EVPA_xy  ,CP_xy   ]
# limits_colors_4panel_xy = [(0,amax(I_xy)),(0,8e-1),(-90.,90.),(0,1e-1)]
# figure(8):              [I_xy    ,|Q,U|_xy,EVPA_xy  ,|V|_xy   ]
# limits_colors_4panel_xy = [(0,amax(I_xy)),(0,amax(abs(Q_xy+1j*U_xy))),(-90.,90.),(0,amax(V_xy))]

# for panel in [0,1,3]:
titles_4panel_xy[0]=titles_4panel_xy[0][:-1]+r" \times 10^{"+str(int(log10(I_xy_scale)))+"}$"
titles_4panel_xy[1]=titles_4panel_xy[1][:-1]+r" \times 10^{"+str(int(log10(LP_xy_scale)))+"}$"
titles_4panel_xy[3]=titles_4panel_xy[3][:-1]+r" \times 10^{"+str(int(log10(V_xy_scale)))+"}$"
titles_miniversion_4panel_uv[0]=titles_miniversion_4panel_uv[0][:-2]+r" \times 10^{"+str(int(log10(I_uv_scale)))+"}|$"

# limits_colors_4panel_xy = [(0,amax(I_xy)),(0,amax(abs(Q_xy+1j*U_xy))),(0.,180.),(amin(V_xy),amax(V_xy))]
limits_colors_4panel_xy = [
    (0,amax(I_xy)*I_xy_scale),
    (0,amax(abs(LP_xy))*LP_xy_scale),
    (0.,180.),
    (amin(V_xy)*V_xy_scale,amax(V_xy)*V_xy_scale)]
# limits_colors_4panel_xy = [(0,round(amax(I_xy)*I_xy_scale,0)),(0,round(amax(abs(Q_xy+1j*U_xy))*I_xy_scale,0)),(0.,180.),(round(amin(V_xy)*I_xy_scale,1),round(amax(V_xy)*I_xy_scale,1))]
# same as above but without round
# limits_colors_4panel_xy = [(0,amax(I_xy)*I_xy_scale),(0,amax(abs(Q_xy+1j*U_xy))*I_xy_scale),(0.,180.),(amin(V_xy)*I_xy_scale,amax(V_xy)*I_xy_scale)]
#limits_colors_4panel_xy = [(0,6e-4),(0,8e-1),(-90.,90.),(0,1e-1)]
# figure(1):       [I_uv    ,Q_uv        ,U_uv        ,V_uv        ]
# limits_colors_uv = [(0,amax(abs(I_uv))*I_uv_scale),(0,amax(abs(Q_uv))),(0,amax(abs(U_uv))),(0,amax(abs(V_uv)))]
limits_colors_uv = [(0,amax(abs(I_uv))),(0,amax(abs(Q_uv))),(0,amax(abs(U_uv))),(0,amax(abs(V_uv)))]
#limits_colors_uv = [(0,3),(0,0.5),(0,0.4),(0,0.1)]
# figure(2): [I_uv,mbreve_uv,EVPA_uv,|V_uv/I_uv|]
limits_colors_miniversion_4panel_uv = [(0,amax(abs(I_uv))*I_uv_scale),(0,1.0),(-90,90),(0,1.)]
# limits_colors_miniversion_4panel_uv = [(0,around(amax(abs(I_uv)),0)),(0,1.0),(-90,90),(0,1.)]
#limits_colors_miniversion_4panel_uv = [(0,3),(0,1.0),(-90,90),(0,0.5)]

limits_xy = array([-image_size/2+20,image_size/2-20,-image_size/2+20,image_size/2-20])

# size_tmp=70;limits_xy = [-size_tmp,size_tmp,-size_tmp,size_tmp] ## FOR 102Ghz case 418

#limits_xy = [-200,200,-200,200] # helical pattern in jet, filaments in QUV
#limits_xy = [-45,45,-45,45]
if observing_frequency>200.:
    limits_uv = [-10,10,-10,10]
else:
    limits_uv = [-6,6,-6,6]

 
# figure(3)
## PF_xy does not asymptote to zero nicely -> windowing
#######################################################
# http://dsp.stackexchange.com/questions/1108/fft-of-image-data-mirroring-to-avoid-boundary-effects
# http://dsp.stackexchange.com/a/362/392
# http://blogs.mathworks.com/steve/2009/12/04/fourier-transform-visualization-using-windowing/
#######################################################
##pcolormesh(u_no_zeropadding,v_no_zeropadding,log10( fftshift(fft2( (Q_xy+1j*U_xy)/I_xy) ,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor])) ) # looks weird, but inverse transform looks good
##pcolormesh(log10( rfft2( abs(Q_xy+1j*U_xy)/I_xy ,s=[int(nxy)*1,int(nxy)*1]) ),cmap=cm.afmhot) # looks weird, but inverse transform looks good
#mtilde_mean=0. # shape=[int(nxy)*4,int(nxy)average(abs(Q_xy+1j*U_xy)/I_xy)
##pcolormesh(u_no_zeropadding,v_no_zeropadding,log10( abs(fftpack.fftshift(fftpack.fft2( abs(Q_xy+1j*U_xy)/I_xy -mtilde_mean ,shape=[int(nxy)*1,int(nxy)*1])) )),cmap=cm.afmhot ) # looks weird, but inverse transform looks good
#pcolormesh(u,v,abs(fftpack.fftshift(fftpack.fft2( abs(Q_xy+1j*U_xy)/I_xy * abs(Q_xy+1j*U_xy) ,shape=[int(nxy)*4,int(nxy)*4])) ) / abs(fftpack.fftshift(fftpack.fft2(abs(Q_xy+1j*U_xy),shape=[int(nxy)*4,int(nxy)*4])) ),cmap=cm.afmhot ) # looks weird, but inverse transform looks good
#plot_shadows("uv")
#colorbar() # format=ticker.FuncFormatter(fmt))
#axis(limits_uv)
#clim(0,10)
#gca().set(title=r"$\tilde{m} \equiv | FFT( (Q+iU)/I ) | \quad (log scale)$",xlabel="$u[G\lambda]$",ylabel="$v[G\lambda]$")
# savefig(filename_out.replace(".png","_mtilde_uv.png"))


################
## IQUV-PLOTS ##
if "IQUV" in WANTED_PLOTS:
    fig_xy = figure(0)
    fig_xy.subplots_adjust(wspace=0.35,hspace=0.22)
    for plot_loop in range(len(titles)):
        
        ## image plane ##
        subplot_inst = fig_xy.add_subplot(221+plot_loop)
        pcolormesh(X,Y,data[:,:,plot_loop]*[I_xy_scale,Q_xy_scale,U_xy_scale,V_xy_scale][plot_loop],cmap=colormaps[plot_loop],norm=[None,MidpointNormalize(midpoint=0)][plot_loop>=1])
        # pcolormesh(X,Y,data[:,:,plot_loop],cmap=colormaps[plot_loop],norm=[None,MidpointNormalize(midpoint=0)][plot_loop>=1])
        cb = colorbar(pad=0)
        # cb = colorbar(format=ticker.FuncFormatter(fmt),pad=0)
        # WIP: I->[0,) Q,U can be negative... mpl does not like (0.0,0,1)
        # cb = colorbar(ticks=(round(limits_colors_xy[plot_loop][0],1),0,round(limits_colors_xy[plot_loop][1],1)),pad=0)
        cb.ax.set_title(r"$Jy/\mu as^2$",fontsize=12)
        # pcolormesh(X,Y,data[:,:,plot_loop],cmap=colormaps[plot_loop],norm=[None,MidpointNormalize(midpoint=0)][plot_loop>=1],vmin=limits_colors_xy[plot_loop][0],vmax=limits_colors_xy[plot_loop][1])
        # colorbar(format=ticker.FuncFormatter(fmt),pad=0,ticks=linspace(limits_colors_xy[plot_loop][0],limits_colors_xy[plot_loop][1],5))
        # #clim(limits_colors_xy[plot_loop])
        
        plot_shadows("xy")
        if plot_loop in [2,3]:
            gca().set(xlabel=r"$x[\mu as]$")
        if plot_loop in [0,2]:
            gca().set(ylabel=r"$y[\mu as]$")
        gca().axis(limits_xy)
        xticks(linspace(round(limits_xy[0],-1),round(limits_xy[1],-1),5))
        yticks(linspace(round(limits_xy[2],-1),round(limits_xy[3],-1),5))
        title(titles[plot_loop]+title_vary_string)

    tight_layout(pad=0.1, w_pad=0.2, h_pad=0) # QUV panel color tick labels need more space due to - sign
    savefig(filename_out.replace(".png","_IQUV_xy.png"))
    
    fig_uv_plane = figure(1)
    fig_uv_plane.subplots_adjust(wspace=0.35,hspace=0.22)

    for plot_loop in range(len(titles)):
        fig_uv_plane.add_subplot(221+plot_loop)
        pcolormesh(u,v,abs([I_uv,Q_uv,U_uv,V_uv][plot_loop]),cmap=colormaps_IQUV_uv[plot_loop],vmin=limits_colors_uv[plot_loop][0],vmax=limits_colors_uv[plot_loop][1])
        colorbar(format=ticker.FuncFormatter(fmt),pad=0,ticks=linspace(limits_colors_uv[plot_loop][0],limits_colors_uv[plot_loop][1],5))
        #clim(limits_colors_uv[plot_loop])

        plot_EHT_uv_tracks()
        plot_shadows("uv")

        if plot_loop in [2,3]:
            gca().set(xlabel=r"$u[G\lambda]$")
        if plot_loop in [0,2]:
            gca().set(ylabel=r"$v[G\lambda]$")
        axis(limits_uv)
        title(titles_IQUV_uv[plot_loop]+title_vary_string)
        #title([r"$\rm\tilde{I}$",r"$\rm\tilde{Q}$",r"$\rm\tilde{U}$",r"$\rm\tilde{V}$"][plot_loop]+title_vary_string)

    tight_layout(pad=0.1, w_pad=0., h_pad=0)
    savefig(filename_out.replace(".png","_IQUV_uv.png"))


    # for plot_loop in range(len(titles)):
    #     fig_uv_plane.add_subplot(221+plot_loop)
    #     pcolormesh(u,v,abs([I_uv,Q_uv,U_uv,V_uv][plot_loop]),cmap=colormaps_IQUV_uv[plot_loop],vmin=limits_colors_uv[plot_loop][0],vmax=limits_colors_uv[plot_loop][1])
    #     colorbar(format=ticker.FuncFormatter(fmt),pad=0,ticks=linspace(limits_colors_uv[plot_loop][0],limits_colors_uv[plot_loop][1],5))
    #     # clim(limits_colors_uv[plot_loop])

    #     plot_EHT_uv_tracks()
    #     plot_shadows("uv")

#####################################################################


if "IP" in WANTED_PLOTS:
    fig_xy_4panel = figure(8)
    fig_xy_4panel.subplots_adjust(wspace=0.35,hspace=0.22)

    for plot_loop in range(len(titles)): #
        
        ## image plane ##
        subplot_inst = fig_xy_4panel.add_subplot(221+plot_loop) # use subplot_inst with quiverkey
    
        # pcolormesh(X,Y,[I_xy,abs(Q_xy+1j*U_xy),EVPA_xy,V_xy][plot_loop])
        pcolormesh(X,Y,[I_xy*I_xy_scale,abs(Q_xy+1j*U_xy)*LP_xy_scale,EVPA_xy,V_xy*V_xy_scale][plot_loop],cmap=colormaps_4panel[plot_loop],norm=[None,MidpointNormalize(midpoint=0)][plot_loop==3],rasterized=True)
        if plot_loop==3:
            # colorbar(ticks=around(arange(round(limits_colors_4panel_xy[plot_loop][0],1),round(limits_colors_4panel_xy[plot_loop][1],1),0.2),1),pad=0)
            colorbar(ticks=(round(limits_colors_4panel_xy[plot_loop][0],1),0,round(limits_colors_4panel_xy[plot_loop][1],1)),pad=0)
        elif plot_loop==2:
            colorbar(ticks=linspace(limits_colors_4panel_xy[plot_loop][0],limits_colors_4panel_xy[plot_loop][1],5),pad=0)
        else:
            cbar_4panel_xy=colorbar(ticks=arange(limits_colors_4panel_xy[plot_loop][0],limits_colors_4panel_xy[plot_loop][1],0.2),pad=0)
            cbar_4panel_xy.locator=matplotlib.ticker.MaxNLocator(5)
            cbar_4panel_xy.ax.set_title(r"$Jy/\mu as^2$",fontsize=12)

        # else:
        #     colorbar(format=ticker.FuncFormatter(fmt),pad=0,ticks=linspace(limits_colors_4panel_xy[plot_loop][0],limits_colors_4panel_xy[plot_loop][1],5))
        clim(limits_colors_4panel_xy[plot_loop])

        if plot_loop==0:
            # plot_polticks()
            # plot_polticks(scale_quiver=3e1*I_img_avg*(observing_frequency/230.)**2.)
            plot_polticks(scale_quiver=3e1*mean(I_xy)*(observing_frequency/230.)**2.) # use subplot_inst somehow
            # plot_polticks(scale_quiver=3e-3) # 4e1*I_img_avg*)
            # plot_polticks(width_quiver=0.01,scale_quiver=20,I_threshold=0.1)
        plot_shadows("xy")
        if plot_loop in [2,3]:
            gca().set(xlabel=r"$x[\mu as]$")
        if plot_loop in [0,2]:
            gca().set(ylabel=r"$y[\mu as]$")
        if plot_loop in [1,3]:
            yticks(visible=False)
        if plot_loop in [0,1]:
            xticks(visible=False)
        axis('equal')
        axis(limits_xy)
        title(titles_4panel_xy[plot_loop]+title_vary_string)

    tight_layout(pad=0, w_pad=0, h_pad=0)
    savefig(filename_out.replace(".png","_I-LP-EVPA-CP_xy.pdf"))


########################################################################



#########################
#manager.window.wm_geometry(fig_pos[plot_loop]) # [WIP: need to understand the arg syntax]
#savefig(filename_out.replace(".png","_"+titles[plot_loop]+".png"))
#################
if "IP" in WANTED_PLOTS and miniversion:

  fig_miniversion_uv_4panel = figure(2)
  fig_miniversion_uv_4panel.subplots_adjust(wspace=0.35,hspace=0.22)

  for plot_loop in range(len(titles)):
    ## uv plane ##
    fig_miniversion_uv_4panel.add_subplot(221+plot_loop)

    # pcolormesh(u,v,abs([I_uv,abs((Q_uv+1j*U_uv)/I_uv),EVPA_uv,abs(V_uv/I_uv)][plot_loop]),cmap=colormaps_4panel[plot_loop])
    # mini-version showing just 4 panels: I, mpl, mbreve, EVPA(visibility version). 
    # RG: Why are ampplitudes of mtilde_uv so large?
    # pcolormesh(u,v,[abs(I_uv),abs(mbreve_uv),EVPA_uv,log10(abs(mtilde_uv))][plot_loop],cmap=colormaps_miniversion_4panel[plot_loop])
    # RG: Replace 4th panel with CP
    pcolormesh(u,v,[abs(I_uv)*I_uv_scale,abs(mbreve_uv),EVPA_uv,abs(V_uv/I_uv)][plot_loop],cmap=colormaps_miniversion_4panel[plot_loop],rasterized=True)

    ## COLORBAR ##
    # if plot_loop==2:
    #     colorbar(ticks=linspace(limits_colors_miniversion_4panel_uv[plot_loop][0],limits_colors_miniversion_4panel_uv[plot_loop][1],5),pad=0)
    if plot_loop==0:
        cbar_miniversion_4panel = colorbar(pad=0,ticks=arange(limits_colors_miniversion_4panel_uv[plot_loop][0],limits_colors_miniversion_4panel_uv[plot_loop][1],0.02))
        cbar_miniversion_4panel.locator = matplotlib.ticker.MaxNLocator(5)
        cbar_miniversion_4panel.update_ticks()
    else:
        colorbar(pad=0,ticks=linspace(limits_colors_miniversion_4panel_uv[plot_loop][0],limits_colors_miniversion_4panel_uv[plot_loop][1],[6,6,5,6][plot_loop]))
        # colorbar(format=ticker.FuncFormatter(fmt),pad=0,ticks=linspace(limits_colors_miniversion_4panel_uv[plot_loop][0],limits_colors_miniversion_4panel_uv[plot_loop][1],5))

    clim(limits_colors_miniversion_4panel_uv[plot_loop])

    plot_EHT_uv_tracks()
    plot_shadows("uv")

    if plot_loop in [2,3]:
        gca().set(xlabel=r"$u[G\lambda]$")
    if plot_loop in [0,2]:
        gca().set(ylabel=r"$v[G\lambda]$")
    if plot_loop in [0,1]:
        xticks(visible=False)
    if plot_loop in [1,3]:
        yticks(visible=False)
    axis('equal')
    axis(limits_uv)
    title(titles_miniversion_4panel_uv[plot_loop]+title_vary_string)

tight_layout(pad=0, w_pad=0, h_pad=0)
savefig(filename_out.replace(".png","_miniversion_4panel_uv.pdf"))



#########################################################################



#[close(window) for window in [0,1,2,5]]



figure(9)

## WIP ##
# I mask?
# edges around ticks (arrow shafts)?
# every=5
# I_threshold=0.1
# scale_quiver=0.2
# width_quiver=0.002

## WORKS WELL FOR 230Ghz (SGR A*)
every=6
I_threshold=0.05 
# scale_quiver=1e-4 # 0.2 * amax(I_xy)  # 0.2  good for LP/I

# scale_quiver=1e-4/I_img_avg # 0.2 * amax(I_xy)  # 0.2  good for LP/I
scale_quiver=3e1/mean(I_xy)*(observing_frequency/230.)**2. # 0.2 * amax(I_xy)  # 0.2  good for LP/I

width_quiver=0.01 # 0.01 good for LP/I
## WORKS WELL FOR 102Ghz (SGR A*)
# every=6
# I_threshold=0.05 
# scale_quiver=1e-2*amax(I_xy)  # 0.2  good for LP/I
# width_quiver=0.01 # 0.01 good for LP/I

# M87 & AVERY's MODEL
# scale_quiver=1e-5

Q_masked=(Q_xy/I_xy)
Q_masked[I_xy < I_threshold*amax(I_xy)] = 0. # None # 0.
U_masked=(U_xy/I_xy)
U_masked[I_xy < I_threshold*amax(I_xy)] = 0. # None # 0.


#figure(10)
#subplot(121)

# quiver_fake_inst= quiver(X[::every],Y[::every],Q_masked[::every,::every],U_masked[::every,::every],headlength=0.,headaxislength=0.,headwidth=0.,color="silver",alpha=0.9,pivot="mid",width=width_quiver/scale_quiver,scale_units="xy",angles="uv",scale=scale_quiver) # ,angles=EVPA_xy[::every,::every])
#figure(10)
NIRVANA=1000
#DEFAULT
#quiver_fake_inst= quiver(X[::every]+NIRVANA,Y[::every]+NIRVANA,(Q_xy/I_xy)[::every,::every],(U_xy/I_xy)[::every,::every],headlength=0.,headaxislength=0.,headwidth=0.,color="silver",alpha=0.9,pivot="mid",width=width_quiver/scale_quiver,scale_units="xy",angles="xy",scale=scale_quiver) # ,angles=EVPA_xy[::every,::every])
quiver_fake_inst= quiver(X[::every]+NIRVANA,Y[::every]+NIRVANA,(Q_xy/I_xy)[::every,::every],(U_xy/I_xy)[::every,::every],headlength=0.,headaxislength=0.,headwidth=0.,color="silver",alpha=0.9,pivot="mid",scale_units="xy",angles="xy") # ,angles=EVPA_xy[::every,::every])


#figure(9)
#cla()
# subplot(121)

# pcolormesh(X,Y,I_xy,cmap=cm.cubehelix)
# colorbar()

# # WORK, but gives weird quiverkey behavior... 
# # Could workaround it by (hidden) dummy plot without angles/angles="uv"
# # I_xy_mask = I_xy>I_threshold*amax(I_xy) # use np.ma() instead...

# #Q_mask=np.ma.less_equal(Q_xy/I_xy,amax(I_xy)*ones(shape(I_xy)))
# # Q_mask=[np.ma.greater_equal(I_xy,I_threshold*amax(I_xy)*ones(shape(I_xy)))][0][::every,::every]
# # U_mask=[np.ma.greater_equal(I_xy,I_threshold*amax(I_xy)*ones(shape(I_xy)))][0][::every,::every]

# #quiver(X[every],Y[every],(Q_xy/I_xy)[Q_mask],(U_xy/I_xy)[U_mask])
# quiver_inst= quiver(X[::every],Y[::every],\
# Q_masked[::every,::every],\
# U_masked[::every,::every],\
# headlength=0.,headaxislength=0.,headwidth=0.,color="silver",alpha=0.9,pivot="mid",width=width_quiver/scale_quiver,scale_units="xy",angles=EVPA_xy[::every,::every],scale=scale_quiver)

# #quiver_inst= quiver(X[::every],Y[::every],(Q_xy/I_xy)[Q_mask],(U_xy/I_xy)[U_mask],headlength=0.,headaxislength=0.,headwidth=0.,color="white",alpha=0.9,pivot="mid",width=0.005,scale_units="xy",angles=EVPA_xy[::every,::every])
# axis((-40,40,-40,40))
# # quiverkey(quiver_fake_inst, 0.3, 0.98, 0.5, r'$50\%$', coordinates='figure', labelpos='W',fontproperties={'weight': 'bold', 'size': 15})

# # STRANGE ERROR MESSAGE WITH CRYPTIC TRACE-BACK
# quiverkey(quiver_fake_inst, 0.7, 0.92, 0.5, r'$50\%$', coordinates='figure', labelpos='W',fontproperties={'weight': 'bold', 'size': 15})
# # quiverkey(quiver_inst, 0.7, 0.92, 0.5, r'$50\%$', coordinates='figure', labelpos='W',fontproperties={'weight': 'bold', 'size': 15})

# gcf().patch.set_facecolor('white')
# # WIP: remove angles ->  solves quiverkey issues...
# quiver(X[::every],Y[::every],(Q_xy/I_xy)[::every,::every],(U_xy/I_xy)[::every,::every],headlength=0.,headaxislength=0.,headwidth=0.,color="white",alpha=0.9,pivot="mid",width=0.005,scale_units="xy") # ,angles=EVPA_xy[::every,::every])
# # WIP: angles="uv" see whether it solves quiverkey issues...


## TESTS ##

# One sanity check is to zero out stokes U to leave only Q. Then the ticks should all be vertical or horizontal. If you make all the Stokes Q's positive, then the ticks should all be vertical.

# Next, zero out all Stokes Q to leave only U. Then the ticks should be 45 degrees East of North or 45 degrees West of North. If you make all the Stokes U's positive, then they should point 45 degrees East of North.

# angles="xy"
# Q_xy[:,:]=0 # Q=0: / and \
# U_xy[:,:]=0 # U=0: | and -
# U_xy[:,:]=0; Q_xy[Q_xy<0]=0. # U=0,Q>=0: |
# Q_xy[:,:]=0; U_xy[U_xy<0]=0. # Q=0,U>=0: \

# angles="uv"
# Q_masked[:,:]=0 # Q=0: vertical
# U_masked[:,:]=0 # U=0: horizontal
# U_masked[:,:]=0; Q_masked[Q_masked<0]=0. # U=0,Q>=0: horizontal
# Q_masked[:,:]=0; U_masked[U_masked<0]=0. # Q=0,U>=0: vertical

# angles=EVPA_xy
# Q_masked[:,:]=0 # Q=0: both vertical & horizontal
# U_masked[:,:]=0 # U=0: all diagonal (both directions "\" and "/" )
# U_masked[:,:]=0; Q_masked[Q_masked<0]=0. # U=0,Q>=0: all "/"
# Q_masked[:,:]=0; U_masked[U_masked<0]=0. # Q=0,U>=0: vertical?




# subplot(122)

## INCONSISTENT WITH ASTRORAY CONVENTIONS
## MICHAEL: This is for EAST=LEFT
# EVPA_xy = pi/2. + 0.5*angle(Q_xy+1j*U_xy)
# x_michael = abs(Q_xy+1j*U_xy)*cos( pi/2 + 0.5*angle(Q_xy+1j*U_xy) )
# y_michael = abs(Q_xy+1j*U_xy)*sin( pi/2 + 0.5*angle(Q_xy+1j*U_xy) )

## CONSISTENT WITH ASTRORAY CONVENTION
## MICHAEL: This is for EAST=RIGHT
EVPA_xy = pi/2. - 0.5*angle(Q_xy+1j*U_xy)
x_michael = abs(Q_xy+1j*U_xy)*cos( pi/2 - 0.5*angle(Q_xy+1j*U_xy) )
y_michael = abs(Q_xy+1j*U_xy)*sin( pi/2 - 0.5*angle(Q_xy+1j*U_xy) )

# x_michael/=I_xy;y_michael/=I_xy

Q_masked=x_michael
Q_masked[I_xy < I_threshold*amax(I_xy)] = None # 0.
U_masked=y_michael
U_masked[I_xy < I_threshold*amax(I_xy)] = None # 0.

# experiment with colormaps and influence on shadow
CMAP_SHADOW=cm.cubehelix
# CMAP_SHADOW=cm.jet ## overemphasizes low data values -> makes shadow appear stronger
# CMAP_SHADOW=cm.gnuplot
# CMAP_SHADOW=cm.gnuplot2
# CMAP_SHADOW=cm.CMRmap # cm.magma #cm.plasma #cm.magma # cm.inferno
# limits_xy = [-100,100,-100,100] # compare to Moscibrodzka+ 2009: RADIATIVE MODELS OF SGR A* FROM GRMHD SIMULATIONS


pcolormesh(X,Y,I_xy*I_xy_scale,cmap=CMAP_SHADOW,vmax=limits_colors_xy[0][1])
# pcolormesh(X,Y,I_xy,cmap=cm.cubehelix,vmax=3e-3)
cb_fig9 = colorbar(pad=0)
cb_fig9.ax.set_title(r"$Jy/\mu as^2$",fontsize=12)

# quiver_inst2 = quiver(X[::every],Y[::every],Q_masked[::every,::every],U_masked[::every,::every],headlength=0.,headaxislength=0.,headwidth=0.,color="silver",alpha=0.5,pivot="mid",width=0.005,angles="xy")

plot_shadows("xy")
# plot_polticks(width_quiver=0.01,scale_quiver=3e-3)
# plot_polticks(width_quiver=0.01,scale_quiver=3e-3*(observing_frequency/230.)**2.)
plot_polticks(width_quiver=0.01,scale_quiver=1e2*mean(LP_xy)*(observing_frequency/230.)**2.)

# strange "fan":
#quiverkey(quiver_inst, 0.7, 0.95, 0.5, r'$50\%$', coordinates='figure', labelpos='W') # ,fontsize=18)# ,fontproperties={'weight': 'bold'})
# WIP
# quiverkey(quiver_inst2, 0.7, 0.92, 0.5, r'$50\%$', coordinates='figure', labelpos='W',fontproperties={'weight': 'bold', 'size': 15})

# quiverkey(quiver_fake_inst, 0.8, 0.98, 0.5, r'$50\%$', coordinates='figure', labelpos='W',fontproperties={'weight': 'bold', 'size': 15})
#quiverkey(quiver_inst, 0.8, 0.98, 0.5, r'$50\%$', coordinates='figure', labelpos='W',fontproperties={'weight': 'bold', 'size': 15})


#quiver(X[::every],Y[::every],mbreve_xy[::every,::every],mbreve_xy[::every,::every],headlength=0.,headaxislength=0.,headwidth=0.,color="white",alpha=0.9,scale_units='xy', angles=EVPA_xy[::every,::every],width=0.005,pivot="mid")

gca().set(xlabel=r"$X[\mu as]$",ylabel=r"$Y[\mu as]$",title=titles_4panel_xy[0]+title_vary_string)
axis('equal')
axis(limits_xy) # +array([+10,-10,+10,-10]))
# tight_layout(pad=0.1, w_pad=0., h_pad=0)
# tight_layout(pad=0, w_pad=0, h_pad=0.95)
tight_layout(rect=(0,0,1,1.02))
savefig(filename_out.replace(".png","_polarization-ticks.png"))



######## DONE ###
print "="*8
print "==DONE=="
print "="*8
######## DONE ###
