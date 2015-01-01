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

pc = scipy.constants.parsec # SI
G = scipy.constants.G # SI
c = scipy.constants.c # SI
Msun = 2e30 # SI
M = 4.3e6 * Msun # SAG A*
rg = G*M/c**2
d_SagA = 8.3e3*pc

rad2microarcsec = 360/(2*pi)*3600*1e6
image_size = 8.15 * (2*rg)/d_SagA * rad2microarcsec
image_size_rad = 8.15 * (2*rg)/d_SagA

## SETTINGS ##
#nxy = 201 # image resolution (nxy+1)x(nxy+1) see ASTRO_RAY_main.cpp
filename = sys.argv[1] # should point to shotimage*.dat file
#filename_out = string.join(filename.split(".")[:-1])+".png" # FIXME: removes dot in filename
filename_out = filename.replace(".dat",".png")

limits = [(0,4e-4),(-1e-4,1e-4),(-1e-4,1e-4),(-5e-5,5e-5)]

try:
    colormap = [cm.gnuplot2,cm.PuOr,cm.bwr,cm.RdBu]
except:
    colormap = [cm.hot,cm.hot,cm.hot,cm.hot]

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

figure(5) # visibility

I_FFT = fftpack.fft2(data[:,:,0])
I_FFT = fftpack.fftshift(I_FFT) # low frequencies in center of image

# FIXME: dim-less image length 1 or 2 ?
# {230.86, 12.2} half-size of the square in a picture plane for each frequency
# 12.2 what units \muarcsec?
pixeldim = image_size/nxy # Specify the linear size of a pixel, in \[Mu]as
X = pixeldim*arange(-round(nxy/2),round(nxy/2)+1)
Y = X[:]

freq_unit=1e-9 # uv plane scale
# Michael:
#                vvvv =2 vvvv
# uvspacing = 1/(pixeldim*nxy*1e-6 / 3600.*pi/180.)
# Jon:
uvspacing = image_size_rad/nxy
FFT_FREQ = fftfreq(int(nxy+1),d=uvspacing)*freq_unit
uv = meshgrid(FFT_FREQ,FFT_FREQ)

I_PSD = np.abs(I_FFT)**2
pcolormesh(unique(FFT_FREQ),unique(FFT_FREQ),log10(abs(I_FFT)))
#imshow(log10(I_PSD))
colorbar()

figure(6) ## Polarization fraction PF = \sqrt{Q^2+U^2}/I ##
pcolormesh(X,Y,sqrt(data[:,:,1]**2+data[:,:,2]**2)/data[:,:,0])
# imshow(sqrt(data[:,:,1]**2+data[:,:,2]**2)/data[:,:,0])
axis((-40,40,-40,40))
xlabel(r"$\mu arcsec$");ylabel(r"$\mu arcsec$")
#xlabel(r"$rad$");ylabel(r"$rad$")
colorbar()
title(r"$\sqrt{Q^2+U^2}/I$")


## PLOT ##
for plot in range(len(titles)):
    figure(plot)
    #imshow(data[:,:,plot],cmap=colormap[plot])
    pcolormesh(X,Y,data[:,:,plot],cmap=colormap[plot])
    colorbar()
    clim(limits[plot])
    xlabel(r"$\mu arcsec$")
    ylabel(r"$\mu arcsec$")
    axis((-40,40,-40,40))
    title(titles[plot])
    manager.window.wm_geometry(fig_pos[plot]) # [WIP: need to understand the arg syntax]
    savefig(filename_out.replace(".png","_"+titles[plot]+".png"))



print "="*8
print "==DONE=="
print "="*8
