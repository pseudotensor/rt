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

## Michaels Notebook ##
# (* Source Parameters *)
# SGRAdec = -(29 + 28.118/3600);

# (* Wavelength (meters) *)
# \[Lambda] = 1.3*10^-3;

# (* Telescope Parameters *)

# SMAvec = {-5464523.400, -2493147.080, 2150611.750}/\[Lambda];
# SMTvec = {-1828796.200, -5054406.800, 3427865.200}/\[Lambda];
# CARMAvec = {-2397431.300, -4482018.900, 3843524.500}/\[Lambda];
# LMTvec = {-768713.9637, -5988541.7982, 2063275.9472}/\[Lambda];
# ALMAvec = {2225037.1851, -5441199.1620, -2479303.4629}/\[Lambda];
# PVvec = {5088967.9000, -301681.6000, 3825015.8000}/\[Lambda];
# PdBIvec = {4523998.40 , 468045.240 , 4460309.760}/\[Lambda];
# SPTvec = {0.0 , 0.0, -6359587.3}/\[Lambda];
# Haystackvec = {1492460, -4457299, 4296835}/\[Lambda];
#######################################################


#import matplotlib 
#matplotlib.use("wx")
import pylab,string,scipy
from pylab import *
from scipy import *
from scipy.constants import *
from scipy import fftpack
import matplotlib.ticker as ticker

angle_unit="arcsec" # "rad"
zeropadding_factor=4 

pc = scipy.constants.parsec # SI
G = scipy.constants.G # SI
c = scipy.constants.c # SI
Msun = 2e30 # SI
M = 4.3e6 * Msun # SAG A*
rg = G*M/c**2
d_SagA = 8.3e3*pc

rad2microarcsec = 360/(2*pi)*3600*1e6
image_size_astroray = 12.2 # @230Ghz
image_size = image_size_astroray * (2*rg)/d_SagA * rad2microarcsec
image_size_rad = image_size_astroray * (2*rg)/d_SagA

shadow_schwarzschild = sqrt(27)*2*rg # diameter
shadow_maximally_spinning = 9./2.*2*rg # diameter


## SETTINGS ##
#nxy = 201 # image resolution (nxy+1)x(nxy+1) see ASTRO_RAY_main.cpp
filename = sys.argv[1] # should point to shotimage*.dat file
#filename_out = string.join(filename.split(".")[:-1])+".png" # FIXME: removes dot in filename
filename_out = filename.replace(".dat",".png")

limits_colors = [(0,4e-4),(-1e-4,1e-4),(-1e-4,1e-4),(-5e-5,5e-5)]
limits_xy = [-50,50,-50,50]
limits_uv = [-10,10,-10,10]

try:
    colormap = [cm.gnuplot2,cm.PuOr,cm.bwr,cm.RdBu]
except:
    colormap = [cm.hot,cm.hot,cm.hot,cm.hot]

def fmt(x, pos):
    '''Customized formatter. Credit goes to the answer in:
    http://stackoverflow.com/questions/25983218/scientific-notation-colorbar-in-matplotlib'''
    a, b = '{:.0e}'.format(x).split('e')
    b = int(b)
    return r'${}\times 10^{{{}}}$'.format(a,b) # how to save the white space between a and \times ?
    #return r'${}\times 10^{{{}}}$'.format(a,b)

## READ-IN ##
fp = open(filename,"rb")
header = fromfile(fp,count=20)
nxy=header[2]
data = fromfile(fp,dtype=float64).reshape(nxy+1,nxy+1,5) # 5 different channels I,Q,U,V
fp.close()

# control figure placement on screen 
manager = get_current_fig_manager()
fig_pos=["+0+0","+500+0","+0+500","+500+500","+250+250"]
titles = ["I","Q","U","V"]
#titles = ["I"]

########################
#figure(5) # visibility #
########################

I_FFT = fftpack.fft2(data[:,:,0],shape=[nxy*zeropadding_factor,nxy*zeropadding_factor])
I_FFT = fftpack.fftshift(I_FFT) # low frequencies in center of image

# find the mag and phase -- shift to put 0 wavenumber at the center
F_mag = abs(I_FFT)
F_phase = angle(I_FFT)
# EVPA_FFT_PHASE = angle(mbreve) + 90./pi ## Michael: EVPA

pixeldim = image_size/nxy # Specify the linear size of a pixel, in \[Mu]as
if angle_unit=="rad":
    pixeldim = image_size_rad/nxy

X = pixeldim*arange(-round(nxy/2),round(nxy/2)+1)
Y = X[:]

freq_unit=1e-9 # uv plane scale
# Michael:
#                vvvv =2 vvvv
# uvspacing = 1/(pixeldim*nxy*1e-6 / 3600.*pi/180.)
# Jon:
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


FFT_FREQ = fftfreq(shape(I_FFT)[0],d=uvspacing)*freq_unit

#pcolormesh(unique(FFT_FREQ),unique(FFT_FREQ),abs(I_FFT)) # ,cmap=cm.hot)
#imshow(log10(I_PSD))
#colorbar()
#clim(0,None)
# axis(limits_uv)
# plot_horizons("uv")
# gca().set(xlabel="$G\lambda$",ylabel="$G\lambda$",title="I")
# savefig(filename_out.replace(".png","_PF_uv.png"))


figure(6) ## Polarization fraction PF = \sqrt{Q^2+U^2}/I ##
PF = sqrt(data[:,:,1]**2+data[:,:,2]**2)/data[:,:,0]
pcolormesh(X,Y,PF)
# imshow(sqrt(data[:,:,1]**2+data[:,:,2]**2)/data[:,:,0])
plot_horizons("xy")
axis(limits_xy)

xlabel(r"$\mu arcsec$");ylabel(r"$\mu arcsec$")
#xlabel(r"$rad$");ylabel(r"$rad$")
colorbar()
title(r"$\sqrt{Q^2+U^2}/I$")
savefig(filename_out.replace(".png","_PF.png"))

figure(2)
#zeropadding_factor=1
PF_uv = fftpack.fft2(PF,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor])
PF_uv = fftpack.fftshift(PF_uv) # low frequencies in center of image
pcolormesh(unique(FFT_FREQ),unique(FFT_FREQ),abs(PF_uv)) # looks weird, but inverse transform looks good
plot_horizons("uv")
colorbar()
axis(limits_uv)
clim(0,10) # WIP: huge outlier somewhere
gca().set(title=r"$\tilde{m}\equiv FFT(\sqrt{Q^2+U^2}/I)$",xlabel="u $(G\lambda)$",ylabel="v $(G\lambda)$")
savefig(filename_out.replace(".png","_PF_uv.png"))

I_xy = data[:,:,0]
Q_xy = data[:,:,1]
U_xy = data[:,:,2]
V_xy = data[:,:,3]

I_uv = fftpack.fftshift(fftpack.fft2(I_xy,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
Q_uv = fftpack.fftshift(fftpack.fft2(Q_xy,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
U_uv = fftpack.fftshift(fftpack.fft2(U_xy,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))
V_uv = fftpack.fftshift(fftpack.fft2(V_xy,shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))

figure(3)
pcolormesh(unique(FFT_FREQ),unique(FFT_FREQ),sqrt(abs(Q_uv)**2+abs(U_uv)**2)/abs(I_uv)) #,cmap=colormap[plot])
colorbar()
plot_horizons("uv")
gca().set(title=r"$\breve{m}\equiv \sqrt{\|FFT(Q)\|^2+\|FFT(U)\|^2}/\|I\|$",xlabel="u $(G\lambda)$",ylabel="v $(G\lambda)$")
axis(limits_uv)
savefig(filename_out.replace(".png","_sqrtQ2U2overI_uv.png"))


figure(4)
pcolormesh(unique(FFT_FREQ),unique(FFT_FREQ),sqrt(abs(Q_uv)**2+abs(U_uv)**2)) #,cmap=colormap[plot])
colorbar()
plot_horizons("uv")
gca().set(title=r"$\breve{\mathrm{P}}\equiv \sqrt{\|FFT(Q)\|^2+\|FFT(U)\|^2}$",xlabel="u $(G\lambda)$",ylabel="v $(G\lambda)$")
axis(limits_uv)
savefig(filename_out.replace(".png","_sqrtQ2U2_uv.png"))

figure(7)
pcolormesh(X,Y,sqrt(abs(Q_xy)**2+abs(U_xy)**2)) #,cmap=colormap[plot])
colorbar(format=ticker.FuncFormatter(fmt))
plot_horizons("xy")
gca().set(title=r"$\mathrm{P}\equiv \sqrt{\|Q\|^2+\|U\|^2}$",xlabel="x $(arcsec)$",ylabel="y $(arcsec)$")
axis(limits_xy)
savefig(filename_out.replace(".png","_sqrtQ2U2_xy.png"))

################
## IQUV-PLOTS ##

fig_image_plane = figure(0)
fig_image_plane.subplots_adjust(wspace = 0.35)
for plot in range(len(titles)):

    ## image plane ##
    fig_image_plane.add_subplot(221+plot)
    #subplot(221+plot)
    pcolormesh(X,Y,data[:,:,plot],cmap=colormap[plot])
    colorbar(format=ticker.FuncFormatter(fmt),pad=0)
    #clim(limits_colors[plot])
    plot_horizons("xy")
    if plot in [2,3]:
        gca().set(xlabel=r"$\mu arcsec$")
    if plot in [0,2]:
        gca().set(ylabel=r"$\mu arcsec$")
    gca().axis(limits_xy)
    title(titles[plot])

savefig(filename_out.replace(".png","_IQUV_image_plane.png"))

fig_uv_plane = figure(1)
fig_uv_plane.subplots_adjust(wspace = 0.35)

for plot in range(len(titles)):
    ## uv plane ##
    fig_uv_plane.add_subplot(221+plot)
    pcolormesh(unique(FFT_FREQ),unique(FFT_FREQ),abs(fftpack.fftshift(fftpack.fft2(data[:,:,plot],shape=[nxy*zeropadding_factor,nxy*zeropadding_factor]))) ) #,cmap=colormap[plot])
    #colorbar()
    colorbar(format=ticker.FuncFormatter(fmt),pad=0)
    #clim(limits_colors[plot])

    plot_horizons("uv")

    if plot in [2,3]:
        gca().set(xlabel=r"u $(G\lambda)$")
    if plot in [0,2]:
        gca().set(ylabel=r"v $(G\lambda)$")
    axis(limits_uv)
    title(titles[plot])

savefig(filename_out.replace(".png","_IQUV_uv_plane.png"))

#manager.window.wm_geometry(fig_pos[plot]) # [WIP: need to understand the arg syntax]
#savefig(filename_out.replace(".png","_"+titles[plot]+".png"))

######## DONE ###
print "="*8
print "==DONE=="
print "="*8
######## DONE ###
