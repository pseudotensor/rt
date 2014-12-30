int solvetemperature (doub rz, const doub zz[],doub ff[],void *pas) {
//solving ODE for electron and proton temperature given "internal energy density" and "density" radial profiles
//ff[0]->Te; ff[1]->Tp
    int indT=ndd-1;
    doub irmin=rtab[0], irmax=1.000001*rtab[indT], 
	     iTsmin=Tstab[0], iTsmax=Tstab[indT];//define minimum and maximum radii and correspondent Ts temperatures=(internal energy density)/density

	//implementation of bisection for radius rz
	doub ra=irmin, rb=irmax, rx;
    int ia=0, ib=indT, ix;
    if((ra<=rz) && (rz<=rb)) {
	    while(ib>ia+1){
		    ix=(ia+ib)>>1;
		    rx=rtab[ix];
		    if(rz<rx) {
				ib=ix;
			} else {
				ia=ix;
			};
	    };
    } else {
        //RG:
        //cout << "irmin,irmax,rx="+ra+irmin+irmax+rx+"\n";
        printf("irmin,irmax,rx=%e,%e,%e\n",irmin,irmax,rx);
        printf("Lookup error 1\n Exiting \n");
        exit(-1);
    };

    //interpolate internal energy density at any radius
	ra=rtab[ia]; rb=rtab[ib];
    doub drman=(rz-ra)/(rb-ra);
    ts[stNt]=(1-drman)*Tstab[ia]+drman*Tstab[ib];
    
	doub Dts=(Tstab[ib]-Tstab[ia])/(rtab[ib]-rtab[ia]),                 //internal energy density gradient
         dcap=2.-0.49/(0.7+kb*zz[0]/me/cc/cc)/(0.7+kb*zz[0]/me/cc/cc),  //heat capacity coefficient dcap=1..2
         rhofun=rhocon*pow(rz/rcut,-rhopo),                             //approximate density at rz
         nuc=8.9e-11*pow(zz[0]/(doub)3e10,(doub)-1.5)*rhofun/1e7,       //Coulomb collisions term. Warning: its normalization depends on BH mass!
		 tail=nuc*2.*mp*rgrav*rgrav*rgrav*rz*rz*rhofun/9e16*(zz[1]-zz[0])+ 2*(sqrt(zz[1])/(sqrt(zz[1])+heat*sqrt(zz[0]))-heat/dcap*sqrt(zz[0])/(sqrt(zz[1])+heat*sqrt(zz[0])))*Dts;
	     //energy flux between electrons and protons
    
	//ODE on electron and proton temperatures
	ff[0]=(2*Dts-tail)/(dcap+1);
    ff[1]=tail+(2*Dts-tail)/(dcap+1);
    return(0);
}
