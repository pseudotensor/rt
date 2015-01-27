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

import pylab,string,scipy
from pylab import *
from scipy import *
from scipy.constants import *
from scipy import fftpack
import matplotlib.ticker as ticker


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

################
## USER SPECS ##

if "vary" in sys.argv: # use that for movies where parameter varies (adds current parameter value to plot titles)
    VARY = string.lower(sys.argv[sys.argv.index("vary")+1])
else:
    VARY="" # can be "","magn_cap","theta","r"

SCATTERING = "ON"
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
image_size_astroray = 12.2 # @230Ghz see sftab array in [ASTRORAY_main.cpp]
image_size = image_size_astroray * (2*rg)/d_SagA * rad2microarcsec
image_size_rad = image_size_astroray * (2*rg)/d_SagA
shadow_schwarzschild = sqrt(27)*2*rg # diameter
shadow_maximally_spinning = 9./2.*2*rg # diameter
#######################################################


## SETTINGS ##
#nxy = 201 # image resolution (nxy+1)x(nxy+1) see ASTRO_RAY_main.cpp
filename = sys.argv[1] # should point to shotimage*.dat file
filename_out = filename.replace(".dat",".png")

## READ-IN ##
fp = open(filename,"rb")
header = fromfile(fp,count=20)
nxy=header[2]
data = fromfile(fp,dtype=float64).reshape(nxy+1,nxy+1,5) 
# 4 different channels I,Q,U,V +1 additional "slot"
# total intensity I
# linearly polarized intensity Q,U
# circularly polarized intensity V
fp.close()

try:
    if VARY=="r":
        rr = 1+( float(filename_out.split(str(int(nxy))+"_")[1].split('.')[0]) -100)/5.
        title_vary_string = " at $r="+str(rr)+"M$"
    elif VARY=="t":
        t = ( float(filename_out.split("_")[0].split('fn')[1])-6100 )*4. *rg/c /60.
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


#limits_colors_xy = [(0,4e-4),(-1e-4,1e-4),(-1e-4,1e-4),(-5e-5,5e-5)] ## IQUV_xy
limits_colors_xy = [(0,2e-4),(-2e-5,2e-5),(-2e-5,2e-5),(-8e-6,8e-6)]
limits_colors_4panel_xy = [(0,4e-5),(0,8e-1),(-90.,90.),(0,1e-1)]
limits_colors_uv = [(0,1),(0,0.1),(0,0.15),(0,0.03)]
limits_colors_4panel_uv = [(0,1.8e-1),(0,1.0),(-90,90),(0,1.0)]

limits_xy = [-50,50,-50,50]
limits_uv = [-10,10,-10,10]

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
    colormaps = [cm.gnuplot2,cm.PuOr,cm.PuOr,cm.PuOr]
    colormaps_4panel = [cm.cubehelix,cm.jet,cm.RdBu_r,cm.hot]
except:
    colormaps = [cm.hot,cm.hot,cm.hot,cm.hot]
    colormaps_4panel = [cm.hot,cm.hot,cm.hot,cm.hot]

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
titles_4panel_xy = [r"$\rm PF\equiv \sqrt{\|Q\|^2+\|U\|^2}$",r"$\rm m_{LP}\equiv \sqrt{\|Q\|^2+\|U\|^2}/I$",r"$\rm EVPA=arctan2(Q,U)\times 90/\pi$",r"$\rm m_{CP}\equiv \|V\|/I$"]
#titles_4panel_xy = ["","","",""]
#titles_4panel_uv = [r"$\tilde{\rm I}$",r"$\breve{\rm m}\equiv \sqrt{\|\tilde{Q}\|^2+\|\tilde{U}\|^2}/\|\tilde{I}\|$",r"$\rm EVPA \equiv phase(\tilde{P}/\tilde{I}) \times 90/\pi$",r"$\breve{\rm CP}\equiv \|\tilde{V}\|/\|\tilde{I}\|$"]
titles_4panel_uv = [r"$\breve{\rm PF}\equiv \|\tilde{Q}+i\tilde{U}\|$",r"$\breve{\rm m}\equiv \|\tilde{Q}+i\tilde{U}\|/\|\tilde{I}\|$",r"$\rm EVPA \equiv phase(\tilde{P}/\tilde{I}) \times 90/\pi$",r"$\breve{\rm CP}\equiv \|\tilde{V}\|/\|\tilde{I}\|$"]

pixeldim = image_size/nxy # Specify the linear size of a pixel, in \[Mu]as
if angle_unit=="rad":
    pixeldim = image_size_rad/nxy

X = pixeldim*arange(-round(nxy/2),round(nxy/2)+1)
#X = pixeldim*arange(-round(nxy/2),round(nxy/2))
Y = X[:]

freq_unit=1e-9 # uv plane scale
uvspacing = image_size_rad/nxy

# 150microarcsec = 1.4 in u-v plane 
# 15 microarcsec = 14  in u-v plane
# for plot labels given the
# nominal scaling of 10^9\lambda for the baseline
uv_schwarzschild = 150.*1.4/(shadow_schwarzschild/2 /d_SagA * rad2microarcsec)
uv_maximally_spinning = 150.*1.4/(shadow_maximally_spinning/2 /d_SagA * rad2microarcsec)

def plot_horizons(domain):
    if string.lower(domain)=="xy":
        gca().add_artist(Circle((0,0),radius = shadow_schwarzschild/2./d_SagA * rad2microarcsec,color="cyan",alpha=0.5,fill=False,lw=4,ls="dashed"))
        gca().add_artist(Circle((0,0),radius = shadow_maximally_spinning/2./d_SagA * rad2microarcsec,color="grey",alpha=0.5,fill=False,lw=4,ls="dashed"))
    elif string.lower(domain)=="uv":
        gca().add_artist(Circle((0,0),radius = uv_schwarzschild,color="cyan",alpha=0.5,fill=False,lw=4,ls="dashed"))
        gca().add_artist(Circle((0,0),radius = uv_maximally_spinning,color="grey",alpha=0.5,fill=False,lw=4,ls="dashed"))
    # return ""



##############
# PROCESSING #
##############

if string.lower(SCATTERING)=="on" or string.lower(SCATTERING)=="yes":
    # r_fresnel = # 300 km [Goodman & Narayan 1989]
    # fresnel_kernel = exp(1j*(r/(2.*r_fresnel))**2.)

    # WIP: Find an appropriate value for 1.3mm SagA*
    # sigma_in_microarcsec = 5.
    # See these references in Fish et al 1409.4690:
    # Bower et al (2006): 1.039 x 0.64 marcsec
    # Shen  et al (2005): 1.41  x 0.75 marcsec
    # @1.3mm point source is blurred to 22muarcsec
    # http://www.cfa.harvard.edu/sma/events/smaConf/posters/images/Doeleman_SMA10_talk.pdf (slide 11 on piercing scattering screen).
    # http://adsabs.harvard.edu/abs/2009A&A...496...77F (Figure 1)
    # http://astro.berkeley.edu/~gbower/ps/sgravlba2.pdf (section on scattering)
    # FWHM(lambda=1cm) = 1mas
    sigma_in_microarcsec = 7.4 # FWHM = 2.3*sigma
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

    for CHANNEL in range(shape(data)[-1]):
        # Or does gaussian_filter only return the kernel?
        data[:,:,CHANNEL] = ndimage.gaussian_filter(data[:,:,CHANNEL],sigma=sigma) ## sigma is in pixel units
    # I_xy = ndimage.gaussian_filter(I_xy,sigma=sigma)
    # Q_xy = ndimage.gaussian_filter(Q_xy,sigma=sigma)
    # U_xy = ndimage.gaussian_filter(U_xy,sigma=sigma)
    # V_xy = ndimage.gaussian_filter(V_xy,sigma=sigma)
    #    CHANNEL *= kernel
    #    CHANNEL *= fresnel_kernel


I_xy = data[:,:,0]
Q_xy = data[:,:,1]
U_xy = data[:,:,2]
V_xy = data[:,:,3]

I_uv = fftpack.fftshift(fftpack.fft2(I_xy,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
Q_uv = fftpack.fftshift(fftpack.fft2(Q_xy,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
U_uv = fftpack.fftshift(fftpack.fft2(U_xy,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
V_uv = fftpack.fftshift(fftpack.fft2(V_xy,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))

u = unique(fftfreq(shape(I_uv)[0],d=uvspacing)*freq_unit)
v = unique(fftfreq(shape(I_uv)[1],d=uvspacing)*freq_unit)

P_uv = fftpack.fftshift(fftpack.fft2(sqrt(Q_xy**2+U_xy**2),shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
Pbreve_xy = abs(Q_xy+1j*U_xy)
Pbreve_uv = abs(Q_uv+1j*U_uv)
#??? mbreve = sqrt(abs(Q_uv)**2+abs(U_uv)**2)/abs(I_uv)
PF_uv = fftpack.fftshift(fftpack.fft2(sqrt(abs(Q_xy)**2+abs(U_xy)**2)/I_xy)) # ,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
#??? mbreve = sqrt(abs(Q_uv)**2+abs(U_uv)**2)/abs(I_uv)

u_no_zeropadding = unique(fftfreq(shape(PF_uv)[0],d=uvspacing)*freq_unit)
v_no_zeropadding = unique(fftfreq(shape(PF_uv)[1],d=uvspacing)*freq_unit)

# find the mag and phase
I_mag = abs(I_uv)
I_phase = angle(I_uv)
EVPA_uv = angle(P_uv/I_uv) * 90./pi ## Michael: EVPA
#EVPA_uv = angle(PF_uv) +*? 90./pi ## Michael: EVPA
#EVPA_uv = angle(mbreve) +*? 90./pi ## Michael: EVPA


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

# figure(9)
# pcolormesh(u,v,abs(P_uv)/abs(I_uv))
# #pcolormesh(u,v,log10(abs(P_uv)/abs(I_uv)),cmap=colormaps[3])
# colorbar(format=ticker.FuncFormatter(fmt))
# clim(-2,2)
# plot_horizons("uv")
# axis(limits_uv)
# gca().set(title=r"$\|\tilde{P}\|/\|\tilde{I}\|$",xlabel="u $(G\lambda)$",ylabel="v $(G\lambda)$")
# savefig(filename_out.replace(".png","_absPtilde_over_absItilde_uv.png"))


# figure(8)
# pcolormesh(u,v,EVPA_uv,cmap=colormaps[3])
# colorbar()
# plot_horizons("uv")
# axis(limits_uv)
# gca().set(title=r"$EVPA \equiv phase(\tilde{P}/\tilde{I})\times 90/\pi$",xlabel="u $(G\lambda)$",ylabel="v $(G\lambda)$")
# savefig(filename_out.replace(".png","_EVPA_uv.png"))


# figure(6) ## Polarization fraction PF = \sqrt{Q^2+U^2}/I ##
# pcolormesh(X,Y,sqrt(Q_xy**2+U_xy**2)/I_xy)
# plot_horizons("xy")
# axis(limits_xy)

# xlabel(r"$\mu arcsec$");ylabel(r"$\mu arcsec$")
# colorbar()
# title(r"$\sqrt{Q^2+U^2}/I$")
# savefig(filename_out.replace(".png","_PF_xy.png"))

# figure(2)
# ## PF_xy data does not asymptote to zero nicely
# ## -> zeropadding no good here, need larger image 
# pcolormesh(u_no_zeropadding,v_no_zeropadding,log10(abs(PF_uv))) # looks weird, but inverse transform looks good
# plot_horizons("uv")
# colorbar(format=ticker.FuncFormatter(fmt))
# axis(limits_uv)
# #clim(0,10) # WIP: huge outlier somewhere
# gca().set(title=r"$\tilde{m}\equiv FFT(\sqrt{Q^2+U^2}/I)$",xlabel="u $(G\lambda)$",ylabel="v $(G\lambda)$")
# savefig(filename_out.replace(".png","_PF_uv.png"))


# figure(3)

# # FPolTab = 
# #   Flatten[ParallelTable[{u, v, 
# #      Abs[(FTelementFast[dataP[[2]], {u*10^9, v*10^9}] + 
# #          I*FTelementFast[dataP[[3]], {u*10^9, v*10^9}])/
# #        FTelementFast[dataP[[1]], {u*10^9, v*10^9}]]}, {u, -10, 10, 
# #      0.5}, {v, -10, 10, 0.5}], 1];

# pcolormesh(u,v,(abs((Q_uv+1j*U_uv)/I_uv)),cmap=colormaps[0])
# titlestring = r"$\breve{m}\equiv \|(\tilde{Q}+i\tilde{U})/\tilde{I}\|$"

# #pcolormesh(u,v,log10(sqrt(abs(Q_uv)**2+abs(U_uv)**2)/abs(I_uv))) #,cmap=colormaps[plot])
# #titlestring = r"$\breve{m}\equiv \sqrt{\|FFT(Q)\|^2+\|FFT(U)\|^2}/\|FFT(I)\|$"
# colorbar(format=ticker.FuncFormatter(fmt))
# clim(0,2)
# plot_horizons("uv")
# gca().set(title=titlestring,xlabel="u $(G\lambda)$",ylabel="v $(G\lambda)$")
# axis(limits_uv)
# savefig(filename_out.replace(".png","_sqrtQ2U2overI_uv.png"))


# figure(4)
# pcolormesh(u,v,sqrt(abs(Q_uv)**2+abs(U_uv)**2)) # ,cmap=colormaps[3])
# colorbar()
# plot_horizons("uv")
# gca().set(title=r"$\breve{\mathrm{P}}\equiv \sqrt{\|\tilde{Q}\|^2+\|\tilde{U}\|^2}$",xlabel="u $(G\lambda)$",ylabel="v $(G\lambda)$")
# axis(limits_uv)
# savefig(filename_out.replace(".png","_sqrtQ2U2_uv.png"))

# figure(7)
# pcolormesh(X,Y,sqrt(abs(Q_xy)**2+abs(U_xy)**2)) #,cmap=colormaps[plot])
# colorbar(format=ticker.FuncFormatter(fmt))
# plot_horizons("xy")
# gca().set(title=r"$\mathrm{P}\equiv \sqrt{\|Q\|^2+\|U\|^2}$",xlabel="x $(\mu arcsec)$",ylabel="y $(\mu arcsec)$")
# axis(limits_xy)
# savefig(filename_out.replace(".png","_sqrtQ2U2_xy.png"))

################
## IQUV-PLOTS ##

fig_xy = figure(0)
fig_xy.subplots_adjust(wspace=0.35,hspace=0.22)
for plot in range(len(titles)):

    ## image plane ##
    fig_xy.add_subplot(221+plot)
    pcolormesh(X,Y,data[:,:,plot],cmap=colormaps[plot],norm=[None,MidpointNormalize(midpoint=0)][plot>=1])
    colorbar(format=ticker.FuncFormatter(fmt),pad=0,ticks=[linspace(limits_colors_xy[plot][0],limits_colors_xy[plot][1],5)])
    clim(limits_colors_xy[plot])
    plot_horizons("xy")
    if plot in [2,3]:
        gca().set(xlabel=r"$x\, \mu arcsec$")
    if plot in [0,2]:
        gca().set(ylabel=r"$y\,\mu arcsec$")
    gca().axis(limits_xy)
    title(titles[plot]+title_vary_string)

savefig(filename_out.replace(".png","_IQUV_xy.png"))

fig_uv_plane = figure(1)
fig_uv_plane.subplots_adjust(wspace=0.35,hspace=0.22)

for plot in range(len(titles)):
    ## uv plane ##
    fig_uv_plane.add_subplot(221+plot)
    pcolormesh(u,v,abs([I_uv,Q_uv,U_uv,V_uv][plot]),cmap=[cm.gnuplot2,cm.jet,cm.jet,cm.jet][plot])
    #pcolormesh(u,v,abs(fftpack.fftshift(fftpack.fft2(data[:,:,plot],shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))) ) #,cmap=colormaps[plot])
    #colorbar()
    colorbar(format=ticker.FuncFormatter(fmt),pad=0,ticks=[linspace(limits_colors_uv[plot][0],limits_colors_uv[plot][1],5)])
    clim(limits_colors_uv[plot])

    plot_horizons("uv")

    if plot in [2,3]:
        gca().set(xlabel=r"u $(G\lambda)$")
    if plot in [0,2]:
        gca().set(ylabel=r"v $(G\lambda)$")
    axis(limits_uv)
    title([r"$\rm\tilde{I}$",r"$\rm\tilde{Q}$",r"$\rm\tilde{U}$",r"$\rm\tilde{V}$"][plot]+title_vary_string)

savefig(filename_out.replace(".png","_IQUV_uv.png"))

#########################

fig_xy_4panel = figure(8)
fig_xy_4panel.subplots_adjust(wspace=0.35,hspace=0.22)

## WIP ##
#? EVPA_xy = angle((Q_xy+1j*U_xy)/I_xy)
EVPA_xy = arctan2(Q_xy,U_xy)*90./pi
#? EVPA_xy = ifftshift(ifft2(EVPA_uv))
#mbreve_xy = sqrt(abs(Q_xy)**2+abs(U_xy)**2)/I_xy
mbreve_xy = abs(Q_xy+1j*U_xy)/abs(I_xy)

#################################
for plot in range(len(titles)): #

    ## image plane ##
    fig_xy_4panel.add_subplot(221+plot)
    pcolormesh(X,Y,[Pbreve_xy,mbreve_xy,EVPA_xy,abs(V_xy)/I_xy][plot],cmap=colormaps_4panel[plot])
    #pcolormesh(X,Y,[data[:,:,plot],abs(Q_xy+1j*U_xy),ifft2(ifftshift(EVPA_uv)),abs(V_xy)][plot],cmap=colormaps[plot])
    if plot==2:
        colorbar(ticks=[linspace(limits_colors_4panel_xy[plot][0],limits_colors_4panel_xy[plot][1],5)])
    else:
        colorbar(format=ticker.FuncFormatter(fmt),pad=0,ticks=[linspace(limits_colors_4panel_xy[plot][0],limits_colors_4panel_xy[plot][1],5)])
    clim(limits_colors_4panel_xy[plot])
    plot_horizons("xy")
    if plot in [2,3]:
        gca().set(xlabel=r"$x\,\mu arcsec$")
    if plot in [0,2]:
        gca().set(ylabel=r"$y\,\mu arcsec$")
    gca().axis(limits_xy)
    title(titles_4panel_xy[plot]+title_vary_string)

savefig(filename_out.replace(".png","_I-LP-EVPA-CP_xy.png"))

#########################

fig_uv_4panel = figure(5)
fig_uv_4panel.subplots_adjust(wspace=0.35,hspace=0.22)

# mbreve_uv = (abs(Q_uv)**2+abs(U_uv)**2)/abs(I_uv)
mbreve_uv = abs(Q_uv+1j*U_uv)/abs(I_uv)
#large_mbreves = sum((mbreve_uv>=1.) * (sqrt(u**2+v**2)<=5.) ) # WIP: within a circle
large_mbreves = sum(mbreve_uv[u<=3.,v<=3.]>=1.) # within a square
print "Fraction of mbreve>1 (within r_uv<3): ",large_mbreves,"/",size(mbreve_uv[u<=3.,v<=3.]),"=",float(large_mbreves) / size(mbreve_uv[u<=3.,v<=3.])

for plot in range(len(titles)):
    ## uv plane ##
    fig_uv_4panel.add_subplot(221+plot)
#    pcolormesh(u,v,abs([I_uv,abs((Q_uv+1j*U_uv)/I_uv),EVPA_uv,abs(V_uv/I_uv)][plot]),cmap=colormaps_4panel[plot])
    pcolormesh(u,v,[Pbreve_uv,mbreve_uv,EVPA_uv,abs(V_uv)/abs(I_uv)][plot],cmap=colormaps_4panel[plot])
    #pcolormesh(u,v,abs(fftpack.fftshift(fftpack.fft2(data[:,:,plot],shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))) ) #,cmap=colormaps[plot])
    #colorbar()
    if plot==2:
        colorbar(ticks=[linspace(limits_colors_4panel_xy[plot][0],limits_colors_4panel_uv[plot][1],5)])
    else:
        colorbar(format=ticker.FuncFormatter(fmt),pad=0,ticks=[linspace(limits_colors_4panel_uv[plot][0],limits_colors_4panel_uv[plot][1],5)])

    clim(limits_colors_4panel_uv[plot])

    plot_horizons("uv")

    if plot in [2,3]:
        gca().set(xlabel=r"u $(G\lambda)$")
    if plot in [0,2]:
        gca().set(ylabel=r"v $(G\lambda)$")
    axis(limits_uv)
    #title("...WIP... r="+str(rr)+"M")
    title(titles_4panel_uv[plot]+title_vary_string)

savefig(filename_out.replace(".png","_4panel_uv.png"))

#########################


#manager.window.wm_geometry(fig_pos[plot]) # [WIP: need to understand the arg syntax]
#savefig(filename_out.replace(".png","_"+titles[plot]+".png"))

######## DONE ###
print "="*8
print "==DONE=="
print "="*8
######## DONE ###
