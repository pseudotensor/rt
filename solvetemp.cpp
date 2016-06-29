int solvetemperature (doub rz, const doub zz[],doub TeTp[],void *pas) {
//solving ODE for electron and proton temperature given "internal energy density" and "density" radial profiles
//TeTp[0]->Te; TeTp[1]->Tp
    int indT=ndd-1;
    doub irmin=rtab[0], irmax=1.000001*rtab[indT], 
	     iTsmin=T_sim_tab[0], iTsmax=T_sim_tab[indT];//define minimum and maximum radii and correspondent Ts temperatures=(internal energy density)/density

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
    ts[stNt]=(1-drman)*T_sim_tab[ia]+drman*T_sim_tab[ib];
    
	doub Dts=(T_sim_tab[ib]-T_sim_tab[ia])/(rtab[ib]-rtab[ia]),                 //internal energy density gradient
         dcap=2.-0.49/(0.7+kb*zz[0]/me/cc/cc)/(0.7+kb*zz[0]/me/cc/cc),  //heat capacity coeTeTpicient dcap=1..2
         rhofun=rhocon*pow(rz/rcut,-rhopo),                             //approximate density at rz
         nuc=8.9e-11*pow(zz[0]/(doub)3e10,(doub)-1.5)*rhofun/1e7,       //Coulomb collisions term. Warning: its normalization depends on BH mass! //RG:FIXME This should be mentinoed in documentation (easy to miss)
		 tail=nuc*2.*mp*rgrav*rgrav*rgrav*rz*rz*rhofun/9e16*(zz[1]-zz[0])+ 2*(sqrt(zz[1])/(sqrt(zz[1])+heat*sqrt(zz[0]))-heat/dcap*sqrt(zz[0])/(sqrt(zz[1])+heat*sqrt(zz[0])))*Dts;
	     //energy flux between electrons and protons
    
	//ODE on electron and proton temperatures
	TeTp[0]=(2*Dts-tail)/(dcap+1);
    TeTp[1]=tail+(2*Dts-tail)/(dcap+1);
    return(0);
}
