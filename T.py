# coding: utf-8
import commands
FILE=sys.argv[-1]
n=array(commands.getoutput("tail -n1 "+FILE).split()[:3],dtype=int)+1
d=loadtxt(FILE,usecols=[3,4,8]).reshape(n[0],n[1],n[2],3)
loglog(d[:,n[1]/2,0,2],mean(d[:,n[1]/2,:,0],axis=1),'c.',label='<T_e>')
loglog(d[:,n[1]/2,0,2],mean(d[:,n[1]/2,:,1],axis=1),'r--',lw=2,label='<T_p>')
legend()
xlabel(r"$r[M]$",fontsize=18)
ylabel(r"$T[K]$",fontsize=18)
savefig("TeTp.png")
