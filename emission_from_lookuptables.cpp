doub lnW=log(nW),             //log of frequency ratio
	 lTe=log(kb*tet/me/cc/cc),//log of electron temperature in rest mass units
	 lTefrac,                 //position in look-up table along temperature axis
	 lTeman,                  //weight of closest look-up table cell along temperature axis
	 lnWfrac,                 //position in look-up table along frequency axis
	 lnWman;                  //weight of closest look-up table cell along frequency axis

if((lTe<lTmin)||(lTe>lTmax)||(lnW<lnWmin)||(lnW>lnWmax)){ // if OUTSIDE of range for emissivities look-up table

	jVc_lookuptab=0.;
	jQc_lookuptab=0.;
	jIc_lookuptab=0.;

 } else { // if INSIDE range of look-up tables use them

    // Te=Tmin -> zeroth element, Te=Tmax -> last element 
	lTefrac=(lTe-lTmin)/(lTmax-lTmin)*Tlen;
	j=floor(lTefrac);
	lTeman=lTefrac-j;

	lnWfrac=(lnW-lnWmin)/(lnWmax-lnWmin)*nWlen;
	k=floor(lnWfrac);
	lnWman=lnWfrac-k;

	temp = lnWman*lTeman;

    //RG: POSSIBLE SPEED-UP by LOCAL VARIABLE FOR EXP() terms
	//look up dimensionless emissivities by 2D multi-linear interpolation
	jIc_lookuptab=exp((-lnWman - lTeman + 1 + temp)*jI[j][k] + (lnWman - temp)*jI[j][1 + k] + (lTeman - temp)*jI[1 + j][k] + temp*jI[1 + j][1 + k]);
	jQc_lookuptab=exp((-lnWman - lTeman + 1 + temp)*jQ[j][k] + (lnWman - temp)*jQ[j][1 + k] + (lTeman - temp)*jQ[1 + j][k] + temp*jQ[1 + j][1 + k]);
	jVc_lookuptab=exp((-lnWman - lTeman + 1 + temp)*jV[j][k] + (lnWman - temp)*jV[j][1 + k] + (lTeman - temp)*jV[1 + j][k] + temp*jV[1 + j][1 + k]);
};

if((lTe<lTminr)||(lTe>lTmax)){                             //if outside of range for rotativities look-up table
	rQc_lookuptab=0.;
	rVc_lookuptab=0.;
} else {
	lTefrac=(lTe-lTminr)/(lTmax-lTminr)*2*Tlen;
	j=floor(lTefrac);
	lTeman=lTefrac-j;

	// look up dimensionless rotativities by 1D interpolation
	rQc_lookuptab=exp((-lTeman + 1)*rQ[j] + (lTeman)*rQ[1 + j]);
	rVc_lookuptab=exp((-lTeman + 1)*rV[j] + (lTeman)*rV[1 + j]);

};





// AGAIN FOR NON-THERMAL

doub lnW_nth=log(nW),             //log of frequency ratio // RG: no nW_nth
//lTe=log(kb*tet/me/cc/cc),//log of electron temperature in rest mass units
	 lTefrac_nth,                 //position in look-up table along temperature axis
	 lTeman_nth,                  //weight of closest look-up table cell along temperature axis
	 lnWfrac_nth,                 //position in look-up table along frequency axis
	 lnWman_nth;                  //weight of closest look-up table cell along frequency axis

if((lTe<lTmin_nth)||(lTe>lTmax_nth)||(lnW_nth<lnWmin_nth)||(lnW_nth>lnWmax_nth)){ // if OUTSIDE of range for emissivities look-up table
  jVc_nth_lookuptab=0.;
  jQc_nth_lookuptab=0.;
  jIc_nth_lookuptab=0.;

 } else { // if INSIDE range of look-up tables use them

  // Te=Tmin -> zeroth element, Te=Tmax -> last element 
  lTefrac_nth=(lTe-lTmin_nth)/(lTmax_nth-lTmin_nth)*Tlen_nth;
  j=floor(lTefrac_nth);
  lTeman_nth=lTefrac_nth-j;

  lnWfrac_nth=(lnW_nth-lnWmin_nth)/(lnWmax_nth-lnWmin_nth)*nWlen;
  k=floor(lnWfrac_nth);
  lnWman_nth=lnWfrac_nth-k;

  temp = lnWman_nth*lTeman_nth;

  //look up dimensionless emissivities by 2D multi-linear interpolation
  jIc_nth_lookuptab=exp((-lnWman_nth - lTeman_nth + 1 + temp)*jI_nth[j][k] + (lnWman_nth - temp)*jI_nth[j][1 + k] + (lTeman_nth - temp)*jI_nth[1 + j][k] + temp*jI_nth[1 + j][1 + k]);
  jQc_nth_lookuptab=exp((-lnWman_nth - lTeman_nth + 1 + temp)*jQ_nth[j][k] + (lnWman_nth - temp)*jQ_nth[j][1 + k] + (lTeman_nth - temp)*jQ_nth[1 + j][k] + temp*jQ_nth[1 + j][1 + k]);
  jVc_nth_lookuptab=exp((-lnWman_nth - lTeman_nth + 1 + temp)*jV_nth[j][k] + (lnWman_nth - temp)*jV_nth[j][1 + k] + (lTeman_nth - temp)*jV_nth[1 + j][k] + temp*jV_nth[1 + j][1 + k]);
};

if((lTe<lTminr_nth)||(lTe>lTmax_nth)){                             //if outside of range for rotativities look-up table
	rQc_nth_lookuptab=0.;
	rVc_nth_lookuptab=0.;
} else {
	lTefrac_nth=(lTe-lTminr_nth)/(lTmax_nth-lTminr_nth)*2*Tlen_nth;
	j=floor(lTefrac_nth);
	lTeman_nth=lTefrac_nth-j;

	//look up dimensionless rotativities by 1D interpolation
	// rQc_nth_lookup=exp((-lTeman_nth + 1)*rQ_nth[j] + (lTeman_nth)*rQ_nth[1 + j]);
	// rVc_nth_lookup=exp((-lTeman_nth + 1)*rV_nth[j] + (lTeman_nth)*rV_nth[1 + j]);

    // LINEAR TABLES (negative values in new scheme)
    if (j==0) cout<<YELLOW"[emission_from_lookuptables.cpp]:"RESET"LINEAR ROTATIVITIES"<<endl;
	rQc_nth_lookuptab=(-lTeman_nth + 1)*rQ_nth[j] + (lTeman_nth)*rQ_nth[1 + j];
	rVc_nth_lookuptab=exp((-lTeman_nth + 1)*rV_nth[j] + (lTeman_nth)*rV_nth[1 + j]);

};
