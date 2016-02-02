int get_electron_temperature (doub T_sim, doub magn, doub& tet, doub& tpt) {

  /****************
   * TEMPERATURES *
   * given T_sim  *
   * get Te,Tp    *
   ****************/

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
  tpt=(1-drman)*tp[ia]+drman*tp[ib];  //compute actual proton and electron temperatures
  tet=(1-drman)*te[ia]+drman*te[ib];
  
  
  /*******************************
   * ADJUST ELECTRON TEMPERATURE *
   *******************************/
  
  if (TEMPERATURE_PRESCRIPTION=="sharma") {
  }
  else if (TEMPERATURE_PRESCRIPTION=="sharma+isoth") {
    doub Te_jet=Te_jet_par*me*cc*cc/kb; // SCS:35 SCS+jet:10
    tet = tet*exp(-magn/magn_cap) + Te_jet*(1.-exp(-magn/magn_cap));
  }
  else if (TEMPERATURE_PRESCRIPTION=="constant_tetp_fraction") {
    tet = ts[ia]/Te_jet_par; // tp/te=3; (assumes ts[ia]=ts[ib])
  }
  else {
    printf(YELLOW"[evalpointzero.cpp]: "RED"NEED TO CHOOSE VALID ELECTRON-TEMPERATURE PRESCRIPTION"RESET"\n");
    exit(1);
  }

  // get_electron_temperature fct
  return(0);
}
