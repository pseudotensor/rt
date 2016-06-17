import numpy,scipy,sys,glob
# import grtrans_rg_rotate.py
from scipy import constants
from numpy import *

def read_astroray_images(filename,fmin=0,fmax=0,view=0):
    ''' READ-IN ASTRORAY IMAGES ...WIP...'''
#    for filename in :
    fp = open(filename,"rb")
    header = numpy.fromfile(fp,count=20)
    fp.close()
    nxy=int(header[2])
    observing_frequency=header[3]
    image_size_inM=header[4]

    pc = scipy.constants.parsec # SI
    G = scipy.constants.G # SI
    c = scipy.constants.c # SI

    Msun = 2e30 # SI
    M = 4.3e6 * Msun # SAG A*
    rg = G*M/c**2
    d_SagA = 8.3e3*pc
    rad2microarcsec = 360/(2*pi)*3600*1e6

    image_size_rad=image_size_inM * (2*rg)/d_SagA
    image_size=image_size_rad * rad2microarcsec
    pixeldim = image_size/nxy # Specify the linear size of a pixel, in \[Mu]as
#    if angle_unit=="rad":
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

    data = zeros((nxy+1,nxy+1,5))

    # COULD REPLACE BY fmin and fmax (see fct args) ...
    # IMAGE_FILES = [entry for entry in sys.argv[1:] if "shotimag" in entry]
    # IMAGE_FILES = glob.glob("astroray-vs-grtrans/shotimag*view*.dat")
    IMAGE_FILES=[filename]

    filename=filename.replace(str(fmin),str(fmin)+"-"+IMAGE_FILES[-1].split("fn")[1].split("case")[0])
    print filename
    for filename_snapshot in IMAGE_FILES:
        fp = open(filename_snapshot,"rb")
        header = fromfile(fp,count=20)
        nxy=int(header[2])
        print "filename: ",filename_snapshot," nxy =",nxy
        data += fromfile(fp,dtype=float64).reshape(nxy+1,nxy+1,5)/size(IMAGE_FILES)
        print "mean(data)=",mean(data)
        fp.close()
    filename_out = filename.replace(".dat",".png")

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

    return data
