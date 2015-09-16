import matplotlib,commands,sys
from commands import getoutput

#matplotlib.use("wx")
try:
    __IPYTHON__ ## ARE WE RUNNING FROM IPYTHON?
except:
    matplotlib.use("Agg") # produces png without X
    from matplotlib import *
    from pylab import *

matplotlib.rc('font', size=20)


## READ-IN ##
F_fastlight=sys.argv[1] # quicka*.dat
F_fdiff20=sys.argv[2]
d_fastlight = loadtxt(F_fastlight)
d_fdiff20 = loadtxt(F_fdiff20)

## PLOT ##
figure(1)
fill_between(d_fastlight[:,0],d_fastlight[:,1],d_fdiff20[:,1],color="red")

## ANNOTATIONS ##
xlabel(r"$\nu/Hz$")
ylabel(r"$F(\nu)/Jy$")
tight_layout()

print "==DONE==\n"
