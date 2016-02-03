# coding: utf-8
d=loadtxt(sys.argv[-1],usecols=[3,4,8]).reshape(272,128,256,3)
loglog(d[:,64,0,2],mean(d[:,64,:,0],axis=1),'k.',label='<T_e>')
loglog(d[:,64,0,2],mean(d[:,64,:,1],axis=1),'r--',lw=2,label='<T_p>')
legend()
xlabel(r"$r[M]$",fontsize=18)
ylabel(r"$T[K]$",fontsize=18)
savefig("TeTp.png")
