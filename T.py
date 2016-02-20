# coding: utf-8
import commands
FILE=sys.argv[-1]
ncut=137 # 120 # 137
rout=3e5 # See Shcherbakov, Penna & McKinney 2012 (section 4.2)
# rout=7.11445236e+03 # WIP
T_fermi=1.5e7 # See Shcherbakov, Penna & McKinney 2012 (section 4.2)

n=array(commands.getoutput("tail -n1 "+FILE).split()[:3],dtype=int)+1
d=loadtxt(FILE,usecols=[3,4,8]).reshape(n[0],n[1],n[2],3)

try:
    # Tsmap=loadtxt("Tsmap-dipole3dfiduaciala.dat") # dipole uses Tsmap.dat
    Tsmap=loadtxt("Tsmap-thermal.dat")
    # Tsmap=loadtxt("Tsmap.dat")
    TSMAP=True
except:
    TSMAP=False
    pass

## PLOTS ##
loglog(d[:,n[1]/2,0,2],mean(d[:,n[1]/2,:,0],axis=1),'c.',label='<T_e>',zorder=1)
loglog(d[:,n[1]/2,0,2],mean(d[:,n[1]/2,:,1],axis=1),'r--',lw=2,label='<T_p>',zorder=1)
if TSMAP:
    loglog(d[:,n[1]/2,0,2],mean(d[:,n[1]/2,:,1],axis=1),'r--',lw=2,label='<T_p>',zorder=1)
scatter(rout,T_fermi,c="b",marker="*",s=100,label="Fermi",zorder=2)
legend(scatterpoints=1,shadow=True,fancybox=True)
axvspan(d[ncut,0,0,2],1e6,fc='grey',alpha=0.25) # annotate radial extension # HARDCODED!
xlabel(r"$r[M]$",fontsize=18)
ylabel(r"$T[K]$",fontsize=18)
tight_layout()
savefig("TeTp.png")
