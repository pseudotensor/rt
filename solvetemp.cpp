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




/*************************************************************************/

int get_electron_temperature (doub T_sim, doub magn, doub ts[], doub te[], doub tp[],doub &tet) {
// int get_electron_temperature (doub T_sim, doub magn, doub ts[], doub te[], doub tp[],doub &tet) {
  // Given one-fluid temperature (as given e.g. by GRMHD simulations)
  // Return proton and electron temperatures (and possibly modify)
  // Assumes solvetemp() (see above in this file) has done its job

  // te=TeTp[0]; // te global variable
  // tp=TeTp[1]; // te global variable
  // doub magn;

  if(T_sim>maxT)                       //if temperature is above allowed, then set it to maximum allowed
	T_sim=maxT;
  if(T_sim<minT)                       //if temperature is below allowed, then set it to minimum allowed
	T_sim=minT;
  int indT=stNt;                      //number of points on temperature look-up grid //RG: WHY INTRODUCE ANOTHER VARIABLE? openMP thread safety?
  doub Ta=ts[0],                      //minimum temperature
    Tb=ts[indT];                   //maximum temperature
  doub Tx,                            //closest temperature on look-up grid
    Tz=T_sim;                       //temperature of interest
  int ia=0;
  int ib=indT;
  if((Ta<=Tz) && (Tz<=Tb)){
	while(ib>ia+1){
      int geod_pt_idx=(ia+ib)>>1;  // RG: bit-shift operations for speed up? (ia+ib)>>1 = (ia+ib)/2 ? golden ratio is faster than 0.5 interval (optimal 1st order method), even better Brent's method
      Tx=ts[geod_pt_idx];
      if(Tz<Tx){
        ib=geod_pt_idx;
      } else {
        ia=geod_pt_idx;
      };
	};
  } else {
	printf(YELLOW"[evalpointzero.cpp]:"RED" Temperature lookup error \nExiting\n"RESET);
	printf(YELLOW"[evalpointzero.cpp]:"RED" Ta=%f Tb=%f Tz=%f T_sim=%f\n"RESET,Ta,Tb,Tz,T_sim);
    //t_solvetrans += (clock() - t_b4_solvetrans) / (doub)CLOCKS_PER_SEC;
	exit(-1);
  };
  Ta=ts[ia];
  Tb=ts[ib];
  doub drman=(Tz-Ta)/(Tb-Ta);              //weight of the closest temperature cell
  doub tpt=(1-drman)*tp[ia]+drman*tp[ib];  //compute actual proton and electron temperatures
  tet=(1-drman)*te[ia]+drman*te[ib];


  /*******************************
   * ADJUST ELECTRON TEMPERATURE *
   *******************************/
  doub Te_jet=Te_jet_par*me*cc*cc/kb; // SCS:35 SCS+jet:10

  // if (magn > magn_cap) 
  //   printf(YELLOW"[evalpointzero.cpp]:"RESET" magn=%e,magn_cap=%e",magn,magn_cap);

  tet = tet*exp(-magn/magn_cap) + Te_jet*(1.-exp(-magn/magn_cap));





  return(0);
}
