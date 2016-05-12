import scipy
from scipy.interpolate import griddata
import commands
from commands import getoutput
rc('font',size=20)

# SPECTRA_FILES="poliresa93.75th*fn5550case0boocr_th.dat" # bad Chi^2: 
SPECTRA_FILES="poliresa93.75th*fn5800case0boocr_th.dat" # min Chi^2: 2.8541
# SPECTRA_FILES="poliresa93.75th*fn*case0boocr_th.dat"

COMBINED_FILE="m_surf_cr_th_781691.dat"
getoutput("grep -h 230.86 "+SPECTRA_FILES+" > "+COMBINED_FILE)
d=loadtxt(COMBINED_FILE)

th=d[:,-1]
x=linspace(amin(th),amax(th),100)
mdot=d[:,-2]
y=linspace(amin(mdot),amax(mdot),100)
X_int,Y_int = meshgrid(x,y)
chisq=d[:,-4]
rhonor=d[:,-5]
C_heat=d[:,-6]

data=d[:,[-1,-2,-4]]

# If your data is like the one you gave in the example, you already have a mesh (you have a value of z for each pair (x,y)) and you only need to reshape the arrays:
# cols = np.unique(xs).shape[0]
# X = xs.reshape(-1, cols)
# Y = ys.reshape(-1, cols)
# Z = zs.reshape(-1, cols)
# answered Jan 25 '14 at 15:11
# Saullo Castro

print "minimum chi^2:",amin(chisq)

CAP=10
chisq[chisq>=CAP]=CAP
chisq[chisq<= 0]= 0

# chisq_grid = griddata(th,mdot,chisq,th,mdot)
# chisq_grid = griddata(th,mdot,chisq,x,y,interp='linear') # mlab
# chisq_grid = griddata((th,mdot),chisq,(x,y),method='linear') # scipy
# chisq_grid = scipy.interpolate.griddata((th,mdot),chisq,(X_int[None,:],Y_int[:,None]),method='linear') # scipy
chisq_grid = scipy.interpolate.griddata((th,mdot),chisq,(X_int,Y_int),method='cubic',fill_value=CAP) # scipy
# scipy.interpolate.griddata((x,z),d,(X_int,Z_int),method='cubic') # scipy
# chisq_grid[chisq_grid>=CAP]=CAP
# chisq_grid[chisq_grid<= 0]= 0

print "minimum gridded chi^2:",amin(chisq_grid)

contours=arange(int(amin(chisq)),11,1)
if size(contours)==0:
    contours=arange(int(amin(chisq)),int(amin(chisq))+11,1)

figure(1)
#C_plot = contour(x/pi*180.,log10(y),chisq_grid,contours,linewidths=3)
C_plot = contour(X_int/pi*180.,Y_int,chisq_grid,contours,linewidths=3)
clabel(C_plot, inline=1, fontsize=15)
# colorbar(pad=0)
xlabel(r"$\Theta / \circ$")
ylabel(r"$\dot{M} / M_\odot/yr$")
tight_layout()

figure(2)
#C_plot = contourf(x/pi*180.,log10(y),chisq_grid,contours)
C_plot = contourf(X_int/pi*180.,Y_int,chisq_grid,contours)
# clabel(C_plot, inline=1, fontsize=12)
colorbar(pad=0)
xlabel(r"$\Theta / \circ$")
ylabel(r"$\dot{M} / M_\odot/yr$")
tight_layout()

getoutput("rm "+COMBINED_FILE)

