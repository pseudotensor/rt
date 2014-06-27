doub lnW=log(nW),             //log of frequency ratio
	 lTe=log(kb*tet/me/cc/cc),//log of electron temperature in rest mass units
	 lTefrac,                 //position in look-up table along temperature axis
	 lTeman,                  //weight of closest look-up table cell along temperature axis
	 lnWfrac,                 //position in look-up table along frequency axis
	 lnWman;                  //weight of closest look-up table cell along frequency axis

if((lTe<lTmin)||(lTe>lTmax)||(lnW<lnWmin)||(lnW>lnWmax)){ //if outside of range for emissivities look-up table
	xVc=0.;
	xQc=0.;
	xIc=0.;
} else {
	lTefrac=(lTe-lTmin)/(lTmax-lTmin)*Tlen;
	j=floor(lTefrac);
	lTeman=lTefrac-j;
	lnWfrac=(lnW-lnWmin)/(lnWmax-lnWmin)*nWlen;
	k=floor(lnWfrac);
	lnWman=lnWfrac-k;
	temp = lnWman*lTeman;
	//look up dimensionless emissivities by 2D multi-linear interpolation
	xIc=exp((-lnWman - lTeman + 1 + temp)*jI[j][k] + (lnWman - temp)*jI[j][1 + k] + (lTeman - temp)*jI[1 + j][k] + temp*jI[1 + j][1 + k]);
	xQc=exp((-lnWman - lTeman + 1 + temp)*jQ[j][k] + (lnWman - temp)*jQ[j][1 + k] + (lTeman - temp)*jQ[1 + j][k] + temp*jQ[1 + j][1 + k]);
	xVc=exp((-lnWman - lTeman + 1 + temp)*jV[j][k] + (lnWman - temp)*jV[j][1 + k] + (lTeman - temp)*jV[1 + j][k] + temp*jV[1 + j][1 + k]);
};

if((lTe<lTminr)||(lTe>lTmax)){                             //if outside of range for rotativities look-up table
	xrQc=0.;
	xrVc=0.;
} else {
	lTefrac=(lTe-lTminr)/(lTmax-lTminr)*2*Tlen;
	j=floor(lTefrac);
	lTeman=lTefrac-j;
	//look up dimensionlesss rotativities by 1D interpolation
	xrQc=exp((-lTeman + 1)*rQ[j] + (lTeman)*rQ[1 + j]);
	xrVc=exp((-lTeman + 1)*rV[j] + (lTeman)*rV[1 + j]);
};