int solvegeodesic(doub t, const doub y[], doub f[], void *params) {
//geodesic solver, code is autogenerated in Wolfram Mathematica (see a correspondent NB file)
//depends on a global variable asq=a*a
//Boyer�Lindquist coordinates are employed

//auxiliary variables
doub rsq=y[4]*y[4], 
	 musq=y[5]*y[5], 
	 sinmu=sqrt(1.-musq),
	 rhosq=rsq + asq*musq, 
	 rhogmin=rsq - asq*musq, 
	 rhosqsq=rhosq*rhosq, 
	 rasq=rsq+asq, 
	 dx=rasq-2.*y[4],
     sinmusq=sinmu*sinmu;

//defining Christoffel symbols=connection coefficients in BL coordinates
doub gg001=(rasq*rhogmin)/(dx*rhosqsq),
	 gg002=(-2.*asq*y[4]*sinmu*y[5])/rhosqsq,
	 gg013=(a*sinmusq*(-3.*rsq*rsq+asq*asq*musq-asq*rsq*(1.+musq)))/(dx*rhosqsq),
	 gg023=(2.*asq*a*y[4]*sinmu*sinmusq*y[5])/rhosqsq,
	 gg100=(dx*rhogmin)/(rhosqsq*rhosq),
	 gg103=-((a*sinmusq*dx*rhogmin)/(rhosqsq*rhosq)),
	 gg111=(1.-y[4])/dx+y[4]/rhosq,
	 gg112=-((asq*sinmu*y[5])/rhosq),
	 gg122=-((y[4]*dx)/rhosq),
	 gg133=-((sinmusq*dx*(rsq*rsq*y[4]+asq*rsq*((doub)-1.+(1.+2.*y[4])*musq)+asq*asq*(musq+((doub)-1.+y[4])*musq*musq)))/(rhosqsq*rhosq)),
	 gg200=(-2.*asq*y[4]*sinmu*y[5])/(rhosqsq*rhosq),
	 gg203=(2.*a*y[4]*rasq*sinmu*y[5])/(rhosqsq*rhosq),
	 gg211=(asq*sinmu*y[5])/(dx*rhosq),
	 gg212=y[4]/rhosq,
	 gg222=-((asq*sinmu*y[5])/rhosq),
	 gg233=-((sinmu*y[5]*(2.*y[4]*rasq*rasq+dx*rhosqsq))/(rhosqsq*rhosq)),
	 gg301=(a*rhogmin)/(dx*rhosqsq),
	 gg302=(-2.*a*y[4]*y[5])/(sinmu*rhosqsq),
	 gg313=((-2.+y[4])*rsq*rsq+asq*rsq*(-1.+(-1.+2.*y[4])*musq)+asq*asq*(musq+(-1.+y[4])*musq*musq))/(dx*rhosqsq),
	 gg323=(y[5]*(2.*asq*y[4]*sinmusq+rhosqsq))/(sinmu*rhosqsq);

//we have a 1-st order ODE on array of size 12
//f[0..7]    : t,  phi,  rp , mup,  r,  mu , tpr, phpr
//y[0..7]=f' : t', phi', r'', mu'', r', mu', t'', ph''

//converting 2-nd order ODE into 1-st order ODE
f[0] = y[6];
f[1] = y[7];
f[4] = y[2];
f[5] = y[3];

//main equations on tpr, rp, mup, phpr
f[6]=-2.*y[2]*y[6]*gg001+(2.*y[3]*y[6]*gg002)/sinmu-2.*y[2]*y[7]*gg013+(2.*y[3]*y[7]*gg023)/sinmu;
f[2]=-(y[6]*y[6]*gg100)-2.*y[6]*y[7]*gg103-y[2]*y[2]*gg111+(2.*y[2]*y[3]*gg112)/sinmu-(y[3]*y[3]*gg122)/(1.-musq)-y[7]*y[7]*gg133;
f[3]=-((y[3]*y[3]*y[5])/(1.-musq))+sinmu*y[6]*y[6]*gg200+2.*sinmu*y[6]*y[7]*gg203+sinmu*y[2]*y[2]*gg211-2.*y[2]*y[3]*gg212+(y[3]*y[3]*gg222)/sinmu+sinmu*y[7]*y[7]*gg233;
f[7]=-2.*y[2]*y[6]*gg301+(2.*y[3]*y[6]*gg302)/sinmu-2.*y[2]*y[7]*gg313+(2.*y[3]*y[7]*gg323)/sinmu;//geodesic equations

//parallel transport of a vector
f[8]=-(y[2]*y[8]*gg001)-y[6]*y[9]*gg001+(y[3]*y[8]*gg002)/sinmu-y[6]*y[10]*gg002-y[7]*y[9]*gg013-y[2]*y[11]*gg013-y[7]*y[10]*gg023+(y[3]*y[11]*gg023)/sinmu;
f[9]=-(y[6]*y[8]*gg100)-y[7]*y[8]*gg103-y[6]*y[11]*gg103-y[2]*y[9]*gg111+(y[3]*y[9]*gg112)/sinmu-y[2]*y[10]*gg112+(y[3]*y[10]*gg122)/sinmu-y[7]*y[11]*gg133;
f[10]=-(y[6]*y[8]*gg200)-y[7]*y[8]*gg203-y[6]*y[11]*gg203-y[2]*y[9]*gg211+(y[3]*y[9]*gg212)/sinmu-y[2]*y[10]*gg212+(y[3]*y[10]*gg222)/sinmu-y[7]*y[11]*gg233;
f[11]=-(y[2]*y[8]*gg301)-y[6]*y[9]*gg301+(y[3]*y[8]*gg302)/sinmu-y[6]*y[10]*gg302-y[7]*y[9]*gg313-y[2]*y[11]*gg313-y[7]*y[10]*gg323+(y[3]*y[11]*gg323)/sinmu;

return(0);
}