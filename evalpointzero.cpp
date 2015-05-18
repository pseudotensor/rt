/**********************************************************************************************/
/* main routine for computing light emission/absorption/propagation properties through plasma */
/**********************************************************************************************/

bool sing;             //whether to do the fast light or full evolution of simulation as the light ray propagates //RG:sing->fast_light
doub rr,               //radius
	 costh,            //cos(polar angle)
	 ph,               //azimuthal angle
	 cossq,            //cos^2(polar angle)
	 sinsq,            //sin^2(polar angle) - to make the code nicer
	 rsq,              //rr*rr
	 Del,              //rsq-2.*rr+asq
	 rhosq,            //rsq+asq*cossq
	 lrman,            //weight of closest cell in interpolation over radial grid
	 thman,            //weight of closest cell in interpolation over theta grid for given radius
	 phman,            //weight of closest cell in interpolation over phi grid
	 phfrac,           //mapping of phi coordinate to a (continuous) cell number in phi direction for fluid simulation
	 timfrac,          //mapping of proper time coordinate to a (continuous) fluid simulation snapshot ID
	 timman,           //weight of the fluid simulation snapshot with the closest ID in interpolation over time
	 temp;             //temp
doub rest[11],         //single record in a fluid simulation file
	 rho,              //physical number density
	 Ttot,             //internal energy density
	 tet,              //proton temperature // RG: typo in comment?? should be electron temperature?
	 tpt;              //electron temperature // RG: typo in comment? should be proton temperature?
doub ly[4],            //coordinate vector along the geodesic
	 kBL[4],           //tangent vector along the geodesic in Boyer-Lindquist (BL) coordinates
	 e1BL[4],          //perpendicular vector in Boyer-Lindquist coordinates parallel-transported along geodesic
	 kupKS[4],         //tangent vector in Kerr-Schild (KS) coordinates
	 e1upKS[4],        //perpendicular vector in Kerr-Schild coordinates
	 kxx[4]={0.,0.,0.,0.},  //tangent vector in locally-flat co-moving reference frame (LFCRF)
	 e1xx[4]={0.,0.,0.,0.}, //perpendicular vector in LFCRF
	 e2xx[4]={0.,0.,0.,0.}, //second perpendicular vector in LFCRF
	 uxx[4]={0.,0.,0.,0.},  //velocity vector in LFCRF  - should be {1.,0.,0.,0.}
     iKS[4][4],        //-,+,+,+ signature covariant (upper indices) Kerr-Schild metric
	 iMKS[4][4],       //modified Kerr-Schild (MKS) metric (code coordinates)
	 MKStoKS[4][4],    //(covariant) vector transformation from MKS to KS
	 u[4],             //4-velocity in MKS
	 ulo[4],           //covariant (low indices) 4-velocity in MKS
	 uKS[4],           //4-velocity in KS
	 uloKS[4],         //covariant (low indices) 4-velocity in KS
	 Bi[4],            //3-vector of magnetic field (Bi[0]==0) (NOT physical magnetic field)
	 Bup[4],           //covariant 4-vector of magnetic field in MKS
	 BupKS[4],         //covariant 4-vector of magnetic field in KS
	 sc[4];            //auxiliary for vector normalization
doub B[4]={0.,0.,0.,0.};    //magnetic field in LFCRF (physical magnetic field), NOT equal to Bi[4]
doub yyloKS[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}; //transformation matrix between KS covariant vector and locally flat co-moving reference frame
doub knorm,            //normalization of spatial part of k vector (while k^\mu k_\mu=0) in LFCRF
	 e1norm,           //normalization of e1 vector spatial part in LFCRF
	 ue1,              //scalar product of u.e1
	 ku,               //scalar product of u.k
	 kbpar,            //scalar product of k.B
	 kbperp,           //length of k x B = \sqrt(B^2-kbpar^2)
	 Be1,              //scalar product of e1.B
	 Be2,              //scalar product of e2.B
	 Bn,               //norm of (physical) magnetic field vector = magnetic field strength
	 sin2k,            //sine of double angle of B vector projected onto e1, e2 plane with respect to (-e1) vector - helps in radiative transfer
	 cos2k,            //correspondent cosine
	 Bnu,              //thermal source term in LFCRF
     fr,               //full redshift/Doppler shift
	 nufr,             //frequency in LFCRF taking redshift into account
	 nW,               //effectively a ratio of observed to cyclotron frequency with some convenient factor
	 jIc,aIc,          //thermal emissivity and absorptivity in total intensity
	 jQc,aQc,          //thermal emissivity and absorptivity in linearly polarized intensity (\sqrt(U^2+Q^2))
     jVc,aVc,          //thermal emissivity and absorptivity in circularly polarized intensity
	 jIc_approx,aIc_approx,        //approximation to thermal emissivity and absorptivity in total intensity
	 rVc,rQc,          //Faraday rotation and conversion coefficients
	 jIc_lookuptab,jQc_lookuptab,jVc_lookuptab,      //dimensionless emissivities computed from the look-up tables
	 XX,               //convenient quantity over which the lookup is conducted for rotativities
	 xrVc,xrQc,        //dimensionless rotativities computed from the look-up tables (still approximate)


	 jIc_nth,aIc_nth,          // non-thermal emissivity and absorptivity in total intensity
	 jQc_nth,aQc_nth,          // non-thermal emissivity and absorptivity in linearly polarized intensity (\sqrt(U^2+Q^2))
     jVc_nth,aVc_nth,          // non-thermal emissivity and absorptivity in circularly polarized intensity
	 jIc_nth_lookuptab,jQc_nth_lookuptab,jVc_nth_lookuptab,      //dimensionless emissivities computed from the look-up tables
	 jIc_approx_nth,aIc_approx_nth,        //approximation to non-thermal emissivity and absorptivity in total intensity
	 rVc_nth,rQc_nth,          //Faraday rotation and conversion coefficients
	 jIc_lookuptab_nth,jQc_lookuptab_nth,jVc_lookuptab_nth,      //dimensionless emissivities computed from the look-up tables
	 XX_nth,               //convenient quantity over which the lookup is conducted for rotativities
	 rVc_nth_lookup,rQc_nth_lookup,        //dimensionless rotativities computed from the look-up tables (still approximate)



	 coskb,            //cosine of angle between k and B
	 sinkb;            //sine of angle between k and B (always above 0)
int  m, k, j, n, i,
	 tn,               //closest cell index (for fluid simulation) in theta direction
	 rn,               //closest cell index (for fluid simulation) in radial direction
	 sn;               //closest ID of fluid simulation snapshot

int curr=*(int*)pas;   //thread number passed by address to current thread
if((curr>=nthreads)||(curr<0)){ //check the bounds of thread number
	printf("Error in thread number \n Exiting \n");
	exit(-1);
};
doub nu=ppy[curr].nu;               //define frequency in current thread
int indx=ppy[curr].indx;            //number of points on current geodesic
doub la=ppy[curr].lamx[0],          //minimum affine parameter (zero)
	 lb=ppy[curr].lamx[indx]+1e-14; //maximum affine parameter (~1 for most geodesics hitting the BH, ~2 for most geodesics leaving to infinity)
doub lz=t,                          //affine parameter of interest = proper time
     lx;                            //affine parameter at the closest point on a geodesic
int  ia=0,                          //min index of a point on a geodesic
	 ib=indx,                       //max index of a point on a geodesic
	 ix;                            //index of the closest point on a geodesic

//given the proper time find where on a geodesic we are
if((lz<=lb) && (la<=lz)){
	while(ib>ia+1){
		ix=(ia+ib)>>1;
		lx=ppy[curr].lamx[ix];
		if(lz<lx){
			ib=ix;
		} else {
			ia=ix;
		};
	};
} else {
	printf("Lookup error of closest geodesic point \n Exiting");
	exit(-1);
};

la=ppy[curr].lamx[ia];              //closest points on a geodesic
lb=ppy[curr].lamx[ib];
doub drman=(lz-la)/(lb-la);         //weight of the closest point on a geodesic

for(m=0;m<4;m++)                    //coordinates on a geodesic for a given proper time
	ly[m]=(1-drman)*ppy[curr].cooxx[m][ia]+drman*ppy[curr].cooxx[m][ib];
for(m=0;m<4;m++)                    //tangent vector to a geodesic at a given proper time
	kBL[m]=((1-drman)*ppy[curr].cooxx[m+4][ia]+drman*ppy[curr].cooxx[m+4][ib])/r0;
for(m=0;m<4;m++)                    //perpendicular vector to a geodesic at a given proper time
	e1BL[m]=(1-drman)*ppy[curr].cooxx[m+8][ia]+drman*ppy[curr].cooxx[m+8][ib];

rr=ly[1];                           //current radius
costh=ly[2];                        //current cos(theta)
ph=dphi-ly[3];                      //current phi, taking into account azimuthal offset and "different" orientations of phi on a geodesic and phi in numerical simulations - check for every simulation!
if((fabs(rr)>100000.) || (rr<1.)) { 
  // RG: 
  // (1) So rr<horizon ... No problem (just stop the geodesic!). Shouldn't abort the code
  // (2) BOUNDARIES SHOULD NOT BE HARDCODED
  // (3) GIVE USER HINT WHAT TO DO?
  printf("Radial coordinate on a geodesic rr=%f is out of bounds \nEXITING\n",rr);
  // RG: FIXME: DOES NOT WORK. 
  // if (rr<1.) rr=rmin;
  exit(-1);
};

doub lrz=log(rr);                  //log of radius, since simulation grid is approximately logarithmic
int indr=ndd-1;                    //max index of the array of radial coordinates
doub lrx,                          //log of radius of the closest radial cell of the fluid simulation
	 lra=coord[0][0][0],           //min log of radius
	 lrb=coord[indr][0][0];        //max log of radius
ia=0;                              //min index of radial coordinates array
ib=indr;                           //max index of radial coordinates array

//find the closest radial cells of the fluid simulation/its radial extension
if((lrz<=lrb) && (lra<=lrz)){
	while(ib>ia+1){
		ix=(ia+ib)>>1;
		lrx=coord[ix][0][0];
		if(lrz<lrx){
			ib=ix;
		} else {
			ia=ix;
		};
	};
} else {
  printf("[evalpointzero.cpp]: Requested radius lrz=%lf coord[%d][0][0]=%lf outside of array of radial coordinates\n...EXITING...\n",lrz,ix,coord[ix][0][0]);
	exit(-1);
};
lra=coord[ia][0][0];               //closest log of radius
lrb=coord[ib][0][0];
rn=ia;                             //closest radial cell index
lrman=(lrz-lra)/(lrb-lra);         //weight of closest cell in interpolation over radial grid

//RG:TODO: LOOK INTO THIS AND DIAGNOSE
//find maximum/minimum theta for a given radius on a simulation grid
const doub critan=-((1.-lrman)*theta[rn][0]+lrman*theta[rn+1][0]); 
//if computed theta is too close to the poles, then set it at the maximum theta on a simulation grid - might lead to image artefacts!
if(costh>critan){
	costh=critan;
};
if(costh<-critan){
	costh=-critan;
};

int indth=thlen-1;                 //maximum index of theta coordinates array
doub costhx,                       //cos(theta) for given radial and theta coordinates
	 tha,                          //closest to costh of interest cos(theta) on the grid
	 thb;                          
ia=0;                              //indices of on cos(theta) grid //RG:?
ib=indth;
//find the closest theta cells of the fluid simulation/its radial extension
while(ib>ia+1){
	ix=(ia+ib)>>1;
	costhx=(1.-lrman)*theta[rn][ix]+lrman*theta[rn+1][ix];
	if(costh<costhx){
		ib=ix;
	} else {
		ia=ix;
	};
};
tha=(1.-lrman)*theta[rn][ia]+lrman*theta[rn+1][ia];
thb=(1.-lrman)*theta[rn][ib]+lrman*theta[rn+1][ib];
drman=(costh-tha)/(thb-tha);       //weight of the closest theta cell in fluid simulation
tn=thlen-ib-1;
thman=1-drman;                     //weight defined in a more uniform way = (1 - weight)

phfrac=fmod(ph*phlen/(2*PI)+0.5,phlen); //closest cell in phi direction - PI factor for half-semisphere and 2*PI factor for the full sphere - check for your simulation!
if(phfrac<0)phfrac+=phlen;              //phi is a periodic coordinate
k=floor(phfrac);                        //true closest closest cell in phi direction
phman=phfrac-k;int ak=(1+k)%phlen;//circular boundary condition 

//compute the closest ID of fluid simulation snapshot + decide whether to do fast light or full simultaneous evolution
if(fdiff==0){                           //fast light approximation
	sn=0;
	sing=true;
} else {                                //simultaneous evolution of fluid simulation & geodesics propagation
	timfrac=(t-mintim)/(maxtim-mintim)*2*fdiff;
	sn=floor(timfrac);
	timman=timfrac-sn;
	sing=false;
}
if(t>=maxtim){                          //for too large time offsets use maximum offset fluid simulation snapshot
	sn=0;
	sing=true;
};
if(t<=mintim){                          //for too small time offsets use maximum offset fluid simulation snapshot
	sn=2*fdiff;
	sing=true;
};

//auxiliary quantities
cossq=costh*costh;
sinsq=1-cossq;
rsq=rr*rr;
rhosq=rsq+asq*cossq;
Del=rsq-2.*rr+asq;
temp=2.*rr/rhosq;

//-,+,+,+ signature covariant KS metric matrix
iKS[0][0]=temp-1.;
iKS[0][1]=temp;
iKS[0][2]=0.;
iKS[0][3]=-a*temp*sinsq;
iKS[1][0]=iKS[0][1];
iKS[1][1]=1.+temp;
iKS[1][2]=0.;
iKS[1][3]=-a*(1.+temp)*sinsq;
iKS[2][0]=iKS[0][2];
iKS[2][1]=iKS[1][2];
iKS[2][2]=rhosq;
iKS[2][3]=0.;
iKS[3][0]=iKS[0][3];
iKS[3][1]=iKS[1][3];
iKS[3][2]=iKS[2][3];
iKS[3][3]=sinsq*(rhosq+asq*(1+temp)*sinsq);

//determine physical quantities at a point with given coordinates
if(rr<=rcut){                          //inside the convergence radius do interpolation
	if(sing)                           //fast light approximation - no time axis
		for(m=0;m<11;m++)              //3D space - multi-linear continuous interpolation over 8 points
			rest[m]=(1-lrman)*(1-thman)*(1-phman)*(*uu[sn])[k][tn][rn][m]+lrman*(1-thman)*(1-phman)*(*uu[sn])[k][tn][rn+1][m]+(1-lrman)*thman*(1-phman)*(*uu[sn])[k][tn+1][rn][m]+lrman*thman*(1-phman)*(*uu[sn])[k][tn+1][rn+1][m]+
			(1-lrman)*(1-thman)*phman*(*uu[sn])[ak][tn][rn][m]+lrman*(1-thman)*phman*(*uu[sn])[ak][tn][rn+1][m]+(1-lrman)*thman*phman*(*uu[sn])[ak][tn+1][rn][m]+lrman*thman*phman*(*uu[sn])[ak][tn+1][rn+1][m];
	else                               //proper simultaneous evolution of simulation and light propagation
		for(m=0;m<11;m++)              //4D space - multi-linear continuous interpolation over 16 points
			rest[m]=(1-timman)*((1-lrman)*(1-thman)*(1-phman)*(*uu[sn])[k][tn][rn][m]+lrman*(1-thman)*(1-phman)*(*uu[sn])[k][tn][rn+1][m]+(1-lrman)*thman*(1-phman)*(*uu[sn])[k][tn+1][rn][m]+
			lrman*thman*(1-phman)*(*uu[sn])[k][tn+1][rn+1][m]+(1-lrman)*(1-thman)*phman*(*uu[sn])[ak][tn][rn][m]+lrman*(1-thman)*phman*(*uu[sn])[ak][tn][rn+1][m]+(1-lrman)*thman*phman*(*uu[sn])[ak][tn+1][rn][m]+
			lrman*thman*phman*(*uu[sn])[ak][tn+1][rn+1][m])+timman*((1-lrman)*(1-thman)*(1-phman)*(*uu[sn+1])[k][tn][rn][m]+lrman*(1-thman)*(1-phman)*(*uu[sn+1])[k][tn][rn+1][m]+
			(1-lrman)*thman*(1-phman)*(*uu[sn+1])[k][tn+1][rn][m]+lrman*thman*(1-phman)*(*uu[sn+1])[k][tn+1][rn+1][m]+(1-lrman)*(1-thman)*phman*(*uu[sn+1])[ak][tn][rn][m]+lrman*(1-thman)*phman*(*uu[sn+1])[ak][tn][rn+1][m]+
			(1-lrman)*thman*phman*(*uu[sn+1])[ak][tn+1][rn][m]+lrman*thman*phman*(*uu[sn+1])[ak][tn+1][rn+1][m]);

	rho=rest[0]*rhonor;                //physical density
	Ttot=rest[1]*mp*cc*cc/3/kb/rest[0];//internal energy density
	//testing k values
	//if((40<k)&&(k<80)){rho/=50.;}//test results - the sense of rotation is clockwise for 1.57<th<3.14
	u[0]=rest[4];                      //Lorentz factor
	u[1]=u[0]*rest[5];                 //u^r
	u[2]=u[0]*rest[6];                 //u^\theta
	u[3]=u[0]*rest[7];                 //u^\phi
	Bi[1]=Bnor*rest[8];                //magnetic field 3-vector
	Bi[2]=Bnor*rest[9];
	Bi[3]=Bnor*rest[10];

	for(i=0;i<4;i++)
		for(k=0;k<4;k++){              //interpolated transformation matrix from MKS to KS
			MKStoKS[i][k]=(1-lrman)*(1-thman)*dxdxp[rn][tn][i][k]+lrman*(1-thman)*dxdxp[rn+1][tn][i][k]+(1-lrman)*thman*dxdxp[rn][tn+1][i][k]+lrman*thman*dxdxp[rn+1][tn+1][i][k];
		}
	for(i=0;i<4;i++)
		for(k=0;k<4;k++)
			iMKS[i][k]=0.;             //initialize metric tensor in MKS
	for(i=0;i<4;i++)
		for(k=0;k<4;k++)
			for(j=0;j<4;j++)
				for(m=0;m<4;m++)       //compute metric tensor in MKS
					iMKS[i][k]+=MKStoKS[m][i]*iKS[m][j]*MKStoKS[j][k];
	for(i=1;i<4;i++){
		uKS[i]=0.;
		for(k=0;k<4;k++)
			uKS[i]+=MKStoKS[i][k]*u[k];//compute 4-velocity in KS coordinates
	}
	                                   //compute 0-th component uKS[0] from normalization to achieve higher precision
	doub expr=fabs((uKS[1]*iKS[0][1] + uKS[3]*iKS[0][3])* (uKS[1]*iKS[0][1] + uKS[3]*iKS[0][3]) - iKS[0][0]*(1. + uKS[1]*uKS[1]*iKS[1][1] + 2.*uKS[1]*uKS[3]*iKS[1][3] + uKS[2]*uKS[2]*iKS[2][2] + uKS[3]*uKS[3]*iKS[3][3]));
	uKS[0]=-((uKS[1]*iKS[0][1] + uKS[3]*iKS[0][3] + sqrt(expr))/iKS[0][0]);
	u[0]=uKS[0];                       //update 0-th component of u as well - simple equality is only achievable in modified KS, where only a spatial transformation is performed
	for(m=0;m<4;m++){
		uloKS[m]=0.;                   //initialize covariant 4-velocity in KS
		for(j=0;j<4;j++)
			uloKS[m]+=iKS[m][j]*uKS[j];//compute covariant 4-velocity in KS
	};
	kupKS[0]=kBL[0]+2*rr*kBL[1]/Del;   //convert k vector from BL to KS, since fluid simulations are in KS/MKS
	kupKS[1]=kBL[1];
	kupKS[2]=kBL[2];
	kupKS[3]=kBL[3]+a*kBL[1]/Del;
	e1upKS[0]=e1BL[0]+2*rr*e1BL[1]/Del;//convert e1 vector from BL to KS
	e1upKS[1]=e1BL[1];
	e1upKS[2]=e1BL[2];
	e1upKS[3]=e1BL[3]+a*e1BL[1]/Del;
	ku=0.;
	ue1=0.;
	for(m=0;m<4;m++){                  //compute scalar products k.u and e1.u
		ue1+=uloKS[m]*e1upKS[m];
		ku+=uloKS[m]*kupKS[m];
	};
	for(m=0;m<4;m++)              
		e1upKS[m]-=ue1/ku*kupKS[m];    //update e1 vector to enforce Lorenz gauge (see Shcherbakov & Huang, 2011 paper)

	for(m=0;m<4;m++){
		ulo[m]=0.;
		for(j=0;j<4;j++)
			ulo[m]+=iMKS[m][j]*u[j];   //compute lower-index 4-velocity in MKS
	};
	
	     //transformation matrix between KS contravariant vector and locally flat co-moving reference frame - requires a correction (see below)
	     //this algorithm doesn't require a correction in BL coordinates. Better algorithm may exist for KS coordinates...
	     //this code is largely duplicated in init.cpp
	doub yyKS[4][4]={{uKS[0], uKS[1], uKS[2], uKS[3]},
		             {uKS[0]*uloKS[1], -(uKS[0]*uloKS[0] + uKS[3]*uloKS[3]), 0, uKS[3]*uloKS[1]},
			         {uloKS[2]*uKS[0], uloKS[2]*uKS[1],uloKS[2]*uKS[2]+1, uloKS[2]*uKS[3]},
				     {-uloKS[3]/uloKS[0], 0, 0, 1}}, 
		 //correction vector
	     xprx[4]={uKS[0]*uloKS[1], -(uKS[0]*uloKS[0] + uKS[3]*uloKS[3]), 0, uKS[3]*uloKS[1]},
		 xloprx[4]={0,0,0,0};

		//compute covariant correction vector
	for(j=0;j<4;j++)
		for(n=0;n<4;n++)
			xloprx[j]+=iKS[j][n]*xprx[n];
		//compute transformation matrix between KS covariant vector and locally flat co-moving reference frame
	for(m=0;m<4;m++)
		for(j=0;j<4;j++)
			for(n=0;n<4;n++)
				yyloKS[m][n]+=iKS[j][n]*yyKS[m][j];
		//normalizing transformation matrices
    for(m=0;m<4;m++){
		sc[m]=0.;
		for(j=0;j<4;j++)
			sc[m]+=yyloKS[m][j]*yyKS[m][j];
		sc[m]=sqrt(fabs(sc[m]));
		for(j=0;j<4;j++){
			yyloKS[m][j]/=sc[m];
			yyKS[m][j]/=sc[m];
		};
	};
		//correcting 0-th vector of yyloKS
	for(j=0;j<4;j++)
		yyloKS[0][j]=-yyloKS[0][j];
		//correcting 1-st vector of yyKS
	for(m=0;m<4;m++){
		yyKS[1][m]=xprx[m];
		for(j=0;j<4;j++)
			yyKS[1][m]-=yyKS[3][m]*xloprx[j]*yyKS[3][j];
	};
		//correcting 1-st vector of yyloKS
	for(j=0;j<4;j++){
		yyloKS[1][j]=0.;
		for(n=0;n<4;n++)
			yyloKS[1][j]+=iKS[n][j]*yyKS[1][n];
	};
		//renormalizing 1-st vectors of yyKS and yyloKS
	sc[1]=0.;
	for(j=0;j<4;j++)
		sc[1]+=yyloKS[1][j]*yyKS[1][j];
	sc[1]=sqrt(fabs(sc[1]));
	for(j=0;j<4;j++){
		yyloKS[1][j]/=sc[1];
		yyKS[1][j]/=sc[1];
	};

		//define 4-vector of magnetic field in MKS
	doub blo=(ulo[1]*Bi[1]+ulo[2]*Bi[2]+ulo[3]*Bi[3]);
	Bup[0]=blo;
	for(m=1;m<4;m++)
		Bup[m]=(Bi[m]+blo*u[m])/u[0];
		//define 4-vector of magnetic field in KS
	for(i=0;i<4;i++){
		BupKS[i]=0.;
		for(k=0;k<4;k++)
			BupKS[i]+=MKStoKS[i][k]*Bup[k];
	}
		//computing co-moving magnetic field
	for(m=0;m<4;m++){
		for(j=0;j<4;j++)
			B[m]+=yyloKS[m][j]*BupKS[j];
	};
	m=1;
} else {                 //outside of the convergence radius do radial extension
	uKS[0]=1.;           //calculations partially repeat those for the inner region
	uKS[1]=0.;
	uKS[2]=0.;
	uKS[3]=0.;
	doub expr=fabs((uKS[1]*iKS[0][1] + uKS[3]*iKS[0][3])* (uKS[1]*iKS[0][1] + uKS[3]*iKS[0][3]) - iKS[0][0]*(1. + uKS[1]*uKS[1]*iKS[1][1] + 2.*uKS[1]*uKS[3]*iKS[1][3] + uKS[2]*uKS[2]*iKS[2][2] + uKS[3]*uKS[3]*iKS[3][3]));
	uKS[0]=-((uKS[1]*iKS[0][1] + uKS[3]*iKS[0][3] + sqrt(expr))/iKS[0][0]);
	u[0]=uKS[0];
	for(m=0;m<4;m++){
		uloKS[m]=0.;
		for(j=0;j<4;j++)
			uloKS[m]+=iKS[m][j]*uKS[j];
	};
	kupKS[0]=kBL[0]+2*rr*kBL[1]/Del;
	kupKS[1]=kBL[1];
	kupKS[2]=kBL[2];
	kupKS[3]=kBL[3]+a*kBL[1]/Del;
	e1upKS[0]=e1BL[0]+2*rr*e1BL[1]/Del;
	e1upKS[1]=e1BL[1];
	e1upKS[2]=e1BL[2];
	e1upKS[3]=e1BL[3]+a*e1BL[1]/Del;
	ku=0.;
	ue1=0.;
	for(m=0;m<4;m++){
		ue1+=uloKS[m]*e1upKS[m];
		ku+=uloKS[m]*kupKS[m];
	};
	for(m=0;m<4;m++)
		e1upKS[m]-=ue1/ku*kupKS[m];

	for(m=0;m<4;m++){
		ulo[m]=0.;
		for(j=0;j<4;j++)
			ulo[m]+=iMKS[m][j]*u[j];
	};

	doub yyKS[4][4]={{uKS[0], uKS[1], uKS[2], uKS[3]},
	                 {uKS[0]*uloKS[1], -(uKS[0]*uloKS[0] + uKS[3]*uloKS[3]), 0, uKS[3]*uloKS[1]},
	                 {uloKS[2]*uKS[0], uloKS[2]*uKS[1],uloKS[2]*uKS[2]+1, uloKS[2]*uKS[3]},
	                 {-uloKS[3]/uloKS[0], 0, 0, 1}},
		xprx[4]={uKS[0]*uloKS[1], -(uKS[0]*uloKS[0] + uKS[3]*uloKS[3]), 0, uKS[3]*uloKS[1]},
		xloprx[4]={0,0,0,0};
	for(j=0;j<4;j++)
		for(n=0;n<4;n++)
			xloprx[j]+=iKS[j][n]*xprx[n];
	for(m=0;m<4;m++)
		for(j=0;j<4;j++)
			for(n=0;n<4;n++)
				yyloKS[m][n]+=iKS[j][n]*yyKS[m][j];

	for(m=0;m<4;m++){
		sc[m]=0.;
		for(j=0;j<4;j++)
			sc[m]+=yyloKS[m][j]*yyKS[m][j];
		sc[m]=sqrt(fabs(sc[m]));
		for(j=0;j<4;j++){
			yyloKS[m][j]/=sc[m];
			yyKS[m][j]/=sc[m];
		};
	};
	for(j=0;j<4;j++)
		yyloKS[0][j]=-yyloKS[0][j];
	for(m=0;m<4;m++){
		yyKS[1][m]=xprx[m];
		for(j=0;j<4;j++)
			yyKS[1][m]-=yyKS[3][m]*xloprx[j]*yyKS[3][j];
	};

	for(j=0;j<4;j++){
		yyloKS[1][j]=0.;
		for(n=0;n<4;n++)
			yyloKS[1][j]+=iKS[n][j]*yyKS[1][n];
	};
	sc[1]=0.;
	for(j=0;j<4;j++)
		sc[1]+=yyloKS[1][j]*yyKS[1][j];
	sc[1]=sqrt(fabs(sc[1]));
	for(j=0;j<4;j++){
		yyloKS[1][j]/=sc[1];
		yyKS[1][j]/=sc[1];
	};
	//"else" code is different from this point
	//radial extension of quantities out of r=rcut sphere
	for(m=0;m<5;m++)                  //using the quantites at the sphere r=rcut we find quantities at larger radius by a power-law extension
		rest[m]=(1-thman)*(1-phman)*uext[k][tn][m]+thman*(1-phman)*uext[k][tn+1][m]+(1-thman)*phman*uext[ak][tn][m]+thman*phman*uext[ak][tn+1][m];
	rho=rest[3]*pow(rr/rcut,-rhopo);  //density is extended with power-law index "-rhopo"
	Ttot=rest[4]*pow(rr/rcut,-1.0);   //temperature is extended with power-law index "-1.0"
	Bpo=0.5*(1.+rhopo);               //magnetic field extension slope - NOT synchronized with command line arguments!
	temp=pow(rr/rcut,-Bpo);           //magnetic field is extended with a power-law index "-0.5*(1+rhopo)" to preserve equipartition
	B[1]=rest[0]*temp;
	B[2]=rest[1]*temp;
	B[3]=rest[2]*temp;
}
///////////////////////
//block of modifications of fluid simulation quantities
double magn=(B[1]*B[1]+B[2]*B[2]+B[3]*B[3])/4/PI/mp/rho/cc/cc; //magnetization

int limit_temperature=0;
if (limit_temperature) {
if(1-fabs(costh)<thlimit)            //if thlimit is above 0, then emissivity in the polar region is effectively set to zero
	Ttot=minT;
if(isBcut)                           //set temperature to zero in high magnetization region
	if(((rr<9) && (magn+20./(9-rg)*(rr-rg))>30)||((rr>=9)&& (magn>10.))) // boundary is in accordance with McKinney et al. (2012)
		Ttot=minT;
if(isBred)                           //temperature reduction in high magnetization regions
	Ttot*=exp(-magn/10.);
 }

int zero_out_rho=1;
if (zero_out_rho) {
  if(1-fabs(costh)<thlimit){            //if thlimit is above 0, then emissivity in the polar region is effectively set to zero
	rho=0.;
    //cout << "Zero out rho due to thlimit...\n";
  }
//double magn=(B[1]*B[1]+B[2]*B[2]+B[3]*B[3])/4/PI/mp/rho/cc/cc; //magnetization
if(isBcut)                           //set temperature to zero in high magnetization region
  //if(((rr<9) && (magn+20./(9-rg)*(rr-rg))>30)||((rr>=9)&& (magn>10.))) {// boundary is in accordance with McKinney et al. (2012)
  if(magn>magn_cap) {// full chop case
		rho=0.;
        //cout << "Zero out rho due isBcut...\n";}
  }
 if(isBred){                           //temperature reduction in high magnetization regions


   /*************************************************************/
   // doub Y_e=0.62; doub mubar=0.62; // solar abundances
   // doub T_e=6.; // =3e10K (T_rest_electron=1=me c^2,me=511kev)

   // // single fluid simulation dump files contains: rho, u, -u^t, -T^r_t/(rho u^r), u^t, v^r, v^theta, v^phi, B^r, B^theta, B^phi
   // T_e = T_{e,code}*exp(-magn/magn_cap) + T_{e,choice} (1-exp(-magn/magn_cap));
   /*************************************************************/

   // limit emission / absorption in jet?
   double rho_jet=1.;
   rho = rho*exp(-magn/magn_cap) + rho*rho_jet*(1.-exp(-magn/magn_cap));

   //rho*=exp(-magn/10.);
   //rho*=exp(-magn/magn_cap); // Limit jet  emission => Lightup disk
   // rho*=(1-exp(-magn/magn_cap)); // Limit disk emission => Lightup jet
  }
}

int theta_slices=0;
if (theta_slices) {
  if ( costh>=cos(trace_theta_slice_angle-trace_theta_slice_width/2.) or costh<=cos(trace_theta_slice_angle+trace_theta_slice_width/2.) ) {
      rho=0.;
    }
 }

int r_slices=0;
if (r_slices) {
  if ( rr<=trace_r_slice-trace_r_slice_width/0.5 or rr>=trace_r_slice+trace_r_slice_width/0.5 ) {
      rho=0.;
    }
 }
///////////////////////


if(Ttot>maxT)                       //if temperature is above allowed, then set it to maximum allowed
	Ttot=maxT;
if(Ttot<minT)                       //if temperature is below allowed, then set it to minimum allowed
	Ttot=minT;
int indT=stNt;                      //number of points on temperature look-up grid
doub Ta=ts[0],                      //minimum temperature
	 Tb=ts[indT];                   //maximum temperature
doub Tx,                            //closest temperature on look-up grid
	 Tz=Ttot;                       //temperature of interest
ia=0;
ib=indT;
if((Ta<=Tz) && (Tz<=Tb)){
	while(ib>ia+1){
		ix=(ia+ib)>>1;
		Tx=ts[ix];
		if(Tz<Tx){
			ib=ix;
		} else {
			ia=ix;
		};
	};
} else {
	printf("Temperature lookup error \n Exiting");
	exit(-1);
};
Ta=ts[ia];
Tb=ts[ib];
drman=(Tz-Ta)/(Tb-Ta);              //weight of the closest temperature cell
tpt=(1-drman)*tp[ia]+drman*tp[ib];  //compute actual proton and electron temperatures
tet=(1-drman)*te[ia]+drman*te[ib];


/*************************************************************/
//RG: ADJUST ELECTRON TEMPERATURE
/*************************************************************/

//doub Y_e=0.62; doub mubar=0.62; // solar abundances
//doub T_e=6.; // =3e10K (T_rest_electron=1=me c^2,me=511kev)

doub Te_jet_par = 35.; // make this a global parameter and include it in fitting procedures

//doub Te_jet_par = 20.; // make this a global parameter and include it in fitting procedures
doub Te_jet=Te_jet_par*me*cc*cc/kb; // SCS:35 SCS+jet:10

//printf("Te_jet=%e,tet=%e\n",Te_jet,tet);

// single fluid simulation dump files contains: rho, u, -u^t, -T^r_t/(rho u^r), u^t, v^r, v^theta, v^phi, B^r, B^theta, B^phi
//doub rho_orig = uu[snapshot_iteration_index][phi_index][theta_index][r_index][0]; // grmhd (HARM) code
//doub u_orig = uu[snapshot_iteration_index][phi_index][theta_index][r_index][1]; // grmhd (HARM) code
//doub gamma_orig=4./3.; // 13./9.; // u_orig/P+1.;
//RG:                               vv
//doub u_Te_cap = Y_e*rho_orig/(mubar*mb)*kb*T_e/(gamma_orig-1.); // = n_e kb T_e/(gamma-1) 
// set u_g=u_Te_cap get T_e the desired isothermal (constant) value, e.g. T_e = 3\times 10^{10}K.;
//(*uu[snapshot_iteration_index])[i][theta_index][r_index][0] = u_orig*(1.-exp(-magn/magn_cap)) + u_Te_cap*(exp(-magn/magn_cap));

// RG: MESS WITH ELECTRON TEMPERATURE

// if (magn > magn_cap) {
//   //printf("magn=%e,magn_cap=%e",magn,magn_cap);
//   //printf("minT=%e",minT);
   // tet = minT; // minT:minimum in temperature table. Tmin is misleading; related to TpTe...
//   tet = Te_jet;
//  }
//tet = minT;

tet = tet*exp(-magn/magn_cap) + Te_jet*(1.-exp(-magn/magn_cap));

//if (tet>1e11) printf("tet=%e\n",tet);
/*************************************************************/
/***********************************************************************************/



for(m=0;m<4;m++)
	for(j=0;j<4;j++)                //geodesic tangential vector in LFCRF 
      kxx[m]+=yyloKS[m][j]*kupKS[j]; // RG: xx?
for(m=0;m<4;m++)
	for(j=0;j<4;j++)                //perpendicular vector in LFCRF
		e1xx[m]+=yyloKS[m][j]*e1upKS[j];
for(m=0;m<4;m++)
	for(j=0;j<4;j++)                //velocity vector in LFCRF (should be {1,0,0,0} )
		uxx[m]+=yyloKS[m][j]*uKS[j];

knorm=0.;
e1norm=0.;                          
for(m=1;m<4;m++){                   //find normalizations of spatial parts of "k" and "e1" vectors
	knorm+=kxx[m]*kxx[m];
	e1norm+=e1xx[m]*e1xx[m];
}
knorm=-sqrt(knorm);                 //since geodesic is computed towards the BH, while radiative transfer is solved backwards, then "k" vector should point in the opposite direction
e1norm=sqrt(e1norm);                //"e1" vector still points in the same direction
for(m=1;m<4;m++){
	kxx[m]/=knorm;
	e1xx[m]/=e1norm;
};
e2xx[0]=0.;                         //find "e2" vector as a vector product e2=k x e1. Note that Melrose defines e2 with a different sign!
e2xx[1]=-e1xx[2]*kxx[3]+e1xx[3]*kxx[2];
e2xx[2]=-e1xx[3]*kxx[1] + e1xx[1]*kxx[3];
e2xx[3]=e1xx[2]*kxx[1] - e1xx[1]*kxx[2];

kbpar=(kxx[1]*B[1]+kxx[2]*B[2]+kxx[3]*B[3]);            //k.B scalar product
kbperp=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]-kbpar*kbpar); //norm of the vector product k x B - main parameter for emissivity

fr=kxx[0];                                              //redshift/Doppler shift 
if(fr<0){                                               //plug for negative frequency. Such problem happens rarely and doesn't interfere with spectrum calculation
  fr=-fr; // RG: HMM... MAYBE "BETTER" TO FLOOR THE FREQUENCY...
	printf("[evalpointzero.cpp]: Negative frequency encountered. Increase accuracy or/and precision. May typically still continue with evaluation...\n");
    printf("(kxx[0],kxx[1],kxx[2],kxx[3]) = (%lf, %lf, %lf, %lf)\n",kxx[0],kxx[1],kxx[2],kxx[3]);
    printf("yyloKS[1][0]*kupKS[0]=%lf, yyloKS[1][1]*kupKS[1]=%lf, yyloKS[2][1]*kupKS[2]=%lf, yyloKS[3][1]*kupKS[3]=%lf\n",yyloKS[1][0]*kupKS[0],yyloKS[1][1]*kupKS[1],yyloKS[2][1]*kupKS[2],yyloKS[3][1]*kupKS[3]); // RG: xx?

    //exit(1);
};


nufr=nu*fr;                                             //frequency in LFCRF taking redshift into account
nW=2*me*cc*2*PI*nufr/(3*ee*kbperp);                     //effectively a ratio of observed to cyclotron frequency with some convenient factor

// PLANCK FCT (Rayleigh-Jeans limit)
Bnu=(2*kb*nufr*nufr*tet)/cc/cc;                         //thermal source term in LFCRF

//computation of dimensionless emissivities from the look-up tables
//RG: actually an interpolation
#include "emission_from_lookuptables.cpp"


//physical emissivities [THERMAL]
//RG: I guess "physical" means cgs units & specific to central mass?
//RG: SEE MPK2012 eq (20) 
doub emissivity_coeff = (sqrt(3.)*ee*ee*ee*rho*rgrav/2.)/(4.*PI*cc*cc*me);

// RG: What about jUc? Obtained via I^2=Q^2+U^2+V^2 ?
// Ioldcoef: \nu q^2/(\sqrt{3} c \nu_\Omega = 4.44\times 10^{-33}\nu/\nu_\Omega
// Qoldcoef: Ioldcoef
// Voldcoef: Ioldcoef (4\cot\theta/3)

// NEW PREFACTORS JON
// I: n_0 \nu
// Q: n_0 \nu
// V: n_0 \nu\cot(\theta)
// alphaI: n_0/\nu
// alphaQ: n_0/\nu
// alphaV: (n_0/\nu)\cot\theta
// rhoQ: n_0/\nu
// rhoV: (n_0/\nu)\cot\theta


// ORIGINAL
// jIc = emissivity_coeff * kbperp * jIc_lookuptab;  //this emissivity is already per unit distance; unit distance is M=rgrav/2
// jQc = emissivity_coeff * kbperp * jQc_lookuptab;
// jVc = fljVc*(ee*ee*ee*rho*rgrav/2.)/(sqrt(3.)*PI*cc*cc*me) * kbpar * jVc_lookuptab;//plus sign for the "k" vector defined above

// doub The=kb*tet/me/cc/cc,                                          //rest mass temperature of electron
// 	 nuc=ee*kbperp/me/cc/2./PI,                                    //cyclotron frequency
// 	 nus=2./9.*nuc*The*The,                                        //peak emissivity frequency
// 	 Xe=nu/nus;                                                    //ratio of observed frequency to peak emissivity frequency

// printf("nu=%e,nufr=%e,Xe=%e\n",nu,nufr,Xe); // nu~230e9,nufr~231e9

jIc = emissivity_coeff * kbperp * jIc_lookuptab;  //this emissivity is already per unit distance; unit distance is M=rgrav/2
jQc = emissivity_coeff * kbperp * jQc_lookuptab; // /4.44e-30*Xe;

jVc = fljVc*(ee*ee*ee*rho*rgrav/2.)/(sqrt(3.)*PI*cc*cc*me) * kbpar * jVc_lookuptab;//plus sign for the "k" vector defined above


//cout<<"[evalpointzero.cpp]:"<<emissivity_coeff<<kbperp<<jIc_nth_lookuptab<<jQc_nth_lookuptab<<kbpar<<jVc_nth_lookuptab<<endl;

// jIc_nth = emissivity_coeff * kbperp * jIc_nth_lookuptab;
// jQc_nth = emissivity_coeff * kbperp * jQc_nth_lookuptab;
// jVc_nth = fljVc * (ee*ee*ee*rho*rgrav/2.)/(sqrt(3.)*PI*cc*cc*me) * kbpar * jVc_nth_lookuptab;//plus sign for the "k" vector defined above
//jVc_nth = jVc; // RG:FIXME SETTING V emissivity equal to old value



// physical emissivities [NON-THERMAL]
// RG: CHECK CONSISTENCY BETWEEN FACTORS/COEFFS HERE VS LOOKUP TABLE
// Jon's tables have different coeffs
// var=OLDSTUFF*(newprefactor/oldprefactor)*newtable
// OLDSTUFF=(sqrt(3.)*ee*ee*ee*rho*rgrav/2.)/(4.*PI*cc*cc*me) * kbperp // ???

// where newprefactor is what I listed as my new prefactor,
// oldprefactor is what *I* list as the old prefactor, and OLDSTUFF is
// what's in astroray just before the variable assignment you showed
// me.

/*********************************************************************
jIc=(sqrt(3.)*ee*ee*ee*rho*rgrav/2.)/(4.*PI*cc*cc*me)*kbperp*jIc_lookuptab;
jQc=(sqrt(3.)*ee*ee*ee*rho*rgrav/2.)/(4.*PI*cc*cc*me)*kbperp*jQc_lookuptab;

4/1/2015 newset2.tgz‎ (1 MB‎) Inbox Wednesday, April 01, 2015 9:15 PM
You replied on 4/1/2015 9:30 PM.  

Roman,

I've created new tables so that there are minimal extra factors to
have in astroray.  That's supposed to minimize mistakes in astroray,
but this is predicated upon you isolated the old functions/factors.

Here are the new function prefactors for the table quantities.  This
applies to the thermal or non-thermal tables I would produce.  Note
that I have yet to produce a thermal table, but I will -- that might
help you isolate factors as well by comparing new and old.  Note these
prefactors exist because they are \nu and \theta and n_0 dependent and
otherwise, e.g., one would have to have another dimension tabulating
\theta or \nu or n_0.

I: n_0 \nu
Q: n_0 \nu
V: n_0 \nu\cot(\theta)
alphaI: n_0/\nu
alphaQ: n_0/\nu
alphaV: (n_0/\nu)\cot\theta
rhoQ: n_0/\nu
rhoV: (n_0/\nu)\cot\theta

The new non-thermal tables and mathematica file are attached.

Here are the original factors to look for in astroray, and one would either:

1) remove prefactors if easy -- best for understanding.  Currently
only remove for non-thermal case, keep for thermal case using old
tables.

or:

2) take my prefactor and divide by these original prefactors I also
provide below (need to avoid singularities though) and don't remove
anything from astroray.  But problem if things like \nu are
dimensionless in some other way than I expect, so then my prefactors
would be used wrongly.

Better is #1, to directly isolate expressions and ensure you see
somewhere the exact coefficients I write.  If you don't see the entire
expression, maybe he's absorbed some factors into making things scaled
a certain way, and then ultimately maybe #2 is fine then.  But best to
understand first.

It depends on how localized the prefactors are and whether you can
readily identify them, understand their normalization, etc.  And
recall the table independent parameters are \nu_\Omega (\nu\Omega in
mathematica), dimensionless T_e, gp (gp replaces T_e as parameter for
non-thermal tables).

Here are the old thermal prefactors:

Ioldcoef: \nu q^2/(\sqrt{3} c \nu_\Omega = 4.44\times 10^{-30}\nu/\nu_\Omega
Qoldcoef: Ioldcoef
Voldcoef: Ioldcoef (4\cot\theta/3)

Recall \alpha's were just using Kirchoff's law inside the code.
rhoQoldcoef,rhoVoldcoef are complicated expressions given by
Shcherbakov 2012 eq 23-26 or Huang & Shcherbakov 58-59.  Ensure you
have isolated each one of those.  I'll write this out just so there's no confusion:

XX= 30.7052T_e/\sqrt{\nu_\Omega}
fx= -0.011 exp(-0.0211864XX) + 2.011exp(-0.212766XX^{1.035})-exp(-0.3663*XX^{1.2})\cos(XX/2)
rhoQoldcoef: 0.00375 fx n_o/(\nu \nu^2_\Omega)
gx=(1-0.11\log(1+0.035XX))
rhoVoldcoef: 0.01126 gx n_0 \cos(\theta)/(\nu\nu_\Omega)

In the above, Te is always the dimensionless temperature.

Recall that the rhoQ,rhoV coefficients are certainly only ever for the thermal case (non-thermal could be done in old prefactor-style, but
makes things not as clean, which is why changed the tables).  But,
that doesn't really matter as you want to isolate these old
prefactors.
****************************************************************************/

// jIc_nth = emissivity_coeff * kbperp * log(exp(jIc_nth_lookuptab)/4.44e-30*Xe);
// jQc_nth = emissivity_coeff * kbperp * log(exp(jQc_nth_lookuptab)/4.44e-30*Xe);
// jVc_nth = fljVc * (ee*ee*ee*rho*rgrav/2.)/(sqrt(3.)*PI*cc*cc*me) * kbpar * jVc_nth_lookuptab;//plus sign for the "k" vector defined above

//RG: too bad stNx undefined here...
//if (stNx==0) printf("jIc=%e,jIc_nth=%e,stNx=%d\n",jIc,jIc_nth,stNx);

//doub newprefactor=rho/me/cc/cc * nufr; // nu:nufr?
//doub oldprefactor=nufr * ee*ee / (sqrt(3.)*cc*nW); // nu:nufr?,nW:nu_omega?
//doub oldprefactor=nW; // DOES NOT WORK
//doub oldprefactor=nW * ee*ee / (sqrt(3.)*cc); // DOES NOT WORK nW:nu_omega/nufr ?
// var=OLDSTUFF*(newprefactor/oldprefactor)*newtable

//jIc_nth = emissivity_coeff * (newprefactor/oldprefactor) * jIc_nth_lookuptab;  //this emissivity is already per unit distance; unit distance is M=rgrav/2
//jIc_nth = rho*rgrav/2. * kbperp * jIc_nth_lookuptab;  //this emissivity is already per unit distance; unit distance is M=rgrav/2
//jQc_nth = emissivity_coeff *            jQc_nth_lookuptab;
//jVc_nth = fljVc            * kbpar    * jVc_nth_lookuptab;//plus sign for the "k" vector defined above


Bn=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);                            //magnetic field strength
coskb=kbpar/Bn;                                                    //cos of angle between "k" and "B" vectors
sinkb=sqrt(1-coskb*coskb);                                         //sin of angle between "k" and "B" vectors
XX=tet/The*sqrt(sqrt(2.)*kbperp*1000.*ee/(2.*cc*me*nufr*PI));      //convenient quantity over which the lookup is conducted for rotativities 

// dimension-full FARADAY CONVERSION and FARADAY ROTATION coefficients [THERMAL]
rQc=flrQc*kbperp*kbperp*rgrav/2.*rho*ee*ee*ee*ee/cc/cc/cc/me/me/me/4/PI/PI/nufr/nufr/nufr*xrQc*(2.011*exp(-pow(XX,(doub)1.035)/4.7)-cos(XX/2.)*exp(-pow(XX,(doub)1.2)/2.73)-0.011*exp(-XX/47.2));
rVc=flrVc*ee*ee*ee*rho*rgrav/2.*kbpar/PI/cc/cc/me/me/nufr/nufr*xrVc*(1.-0.11*log(1.+0.035*XX));


// if (nth && [one point only]) printf("[evalpointzero]: ZERO-OUT NON-THERMAL ROTATIVITIES\n");
// rQc_nth = 0.;
// rVc_nth = 0.;
// RG:FIXME TEMPORARILY SET ROTATIVITIES TO OLD THERMAL VALUES
static int onetime=1;if(onetime) cout<<"[WIP]:...SETTING ROTATION COEFFS TO OLD THERMAL..."<<endl;onetime *= 0;
rQc_nth = rQc;
rVc_nth = rVc;



jIc_nth = emissivity_coeff * kbperp * jIc_nth_lookuptab;
jQc_nth = emissivity_coeff * kbperp * jQc_nth_lookuptab;
jVc_nth = fljVc * (ee*ee*ee*rho*rgrav/2.)/(sqrt(3.)*PI*cc*cc*me) * kbpar * jVc_nth_lookuptab;//plus sign for the "k" vector defined above

// jIc_nth = jIc; // when not resetting I,Q and V ->nan why!?
// jQc_nth = jQc;
// jVc_nth = jVc;


//physical absorptivities
//char absorption[16] = "kirchhoff"; // RG: MAKE THIS A USER CHOICE e.g. inside win_lin_Jon.c
const string absorption = "kirchhoff"; // RG: MAKE THIS A USER CHOICE e.g. inside win_lin_Jon.c
//const string absorption = "tables"; // RG: MAKE THIS A USER CHOICE e.g. inside win_lin_Jon.c
if (absorption=="kirchhoff"){ //FIXME
// if (true){
  aIc=jIc/Bnu;
  aQc=jQc/Bnu;
  aVc=jVc/Bnu;

  // AVOID EXTREMELY SMALL VALUES 
  if (jIc_nth < 1e-32) jIc_nth = 0.; //0.;
  if (jQc_nth < 1e-32) jQc_nth = 0.; //0.;
  if (jVc_nth < 1e-32) jVc_nth = 0.; //0.;

  //RG:FIXME PROBLEM HERE...
  aIc_nth=jIc_nth/Bnu;
  aQc_nth=jQc_nth/Bnu;
  aVc_nth=jVc_nth/Bnu;
  
  // RG: aIc_nth,aQc_nth,aVc_nth are finite

  // AVOID EXTREMELY SMALL VALUES
  if (aIc_nth < 1e-32) aIc_nth = 0.;//aIc; //0.;
  if (aQc_nth < 1e-32) aQc_nth = 0.;//aQc; //0.;
  if (aVc_nth < 1e-32) aVc_nth = 0.;//aVc; //0.;

  if (aIc_nth<0) cout<<"aIc_nth="<<aIc_nth<<endl;

  //if ( max(abs(aIc_nth/aIc),abs(aIc/aIc_nth)) > 2. ) printf("absorption: IQV_old=%e %e %e IQV_new=%e %e %e\n",aIc,aQc,aVc,aIc_nth,aQc_nth,aVc_nth);
  //RG:PROBLEM GOES AWAY WHEN...
  aIc_nth=aIc;
  aQc_nth=aQc;
  aVc_nth=aVc;
 }
 else if (absorption=="powerlaw_jon" or absorption=="tables"){ // See Jon's notebook polarization_loop.nb in ASTRORAY repo

   //cout << "absorption=powerlaw_jon unsupported for now :-(..." << endl;

   // powerlaw index in distribution function f~n^-(b+1) -> N~n^(-b)
   double b=2.5; //RG: be DRY used to defined in init.cpp make it global?
   double gamma_min=100.;
   double omega=2.*PI*nu;
   double Omega_0=ee *Bn/me/cc; // cyclotron frequency = eB/m_e c
 }

 else if (absorption=="powerlaw_approximate"){ // Sazonov (1969) expressions, see eq(52) Huang & Shcherbakov 2011
   // general absorption coefficients eqs(7) Sazonov (1969):
   // eta_I = Im(alpha^{22}+alpha^{11})/nu
   // eta_Q = Im(alpha^{11}-alpha^{22})/nu
   // eta_V = 2Re(alpha^{12})/nu
   cout << "absorption=powerlaw_approximate unavailable for now :-(...better die..." << endl;

   // powerlaw index in distribution function f~n^-(b+1) -> N~n^(-b)
   double b=2.5; //RG: be DRY used to defined in init.cpp make it global?
   double gamma_min=100.;
   double omega=2.*PI*nu;
   double Omega_0=ee *Bn/me/cc; // cyclotron frequency = eB/m_e c

   /* WIP: CHECK EQS */
   //aIc=?;
   // WHERE IS THE EQUATION FOR alpha_I in Huang & Shcherbakov (2011) or Sazonov (1969) ?
   // Take approximate expression for alpha_N =? alpha_I from Avery's PhD thesis see eqn (5.47)
   aIc=sqrt(3.)/24./cc /*WIP-> *f*omega_P**2 */ *omega /*<-WIP*/ *pow(3.*Omega_0/omega*sinkb,(b+2)/2.); // see eqn (5.47) Broderick PhD thesis (2004) assume plasma frequency ~ nu
   aQc=3.1e-4*(b+2)*tgamma((3.*b+2)/12.)*tgamma((3.*b+10)/12.)*(b-1.)/pow(gamma_min,1.-b)*pow(3.*Omega_0/omega*sinkb,(b+2)/2.)/nu; // Gamma fct: tgamma() needs cmath
   aVc=4.1e-4*(b+3)/(b+1)*tgamma((3.*b+7)/12.)*tgamma((3.*b+11)/12.)*(b-1.)/pow(gamma_min,1.-b)*coskb/sinkb*pow(3.*Omega_0/omega*sinkb,(b+3)/2.)/nu; // Gamma fct: tgamma() needs cmath;

   exit(1);
 }
else{
  cout << "[evalpointzero.cpp]: CHOOSE YOUR ABSORPTION. Possible choices are: absorption=\"kirchhoff\" or absorption=\"powerlaw_approximate\" (yet unsupported)"<<endl;
   exit(1);
}


Be1=B[1]*e1xx[1]+B[2]*e1xx[2]+B[3]*e1xx[3];                        //scalar product of B.e1
Be2=B[1]*e2xx[1]+B[2]*e2xx[2]+B[3]*e2xx[3];                        //scalar product of B.e2
sin2k=-2.*Be1*Be2/(Be1*Be1+Be2*Be2);                               //sine of double angle of B vector projected onto e1, e2 plane with respect to (-e1) vector - helps in radiative transfer
cos2k=(Be2*Be2-Be1*Be1)/(Be1*Be1+Be2*Be2);                         //cosine of double angle of B vector projected onto e1, e2 plane with respect to (-e1) vector - helps in radiative transfer

doub The=kb*tet/me/cc/cc,                                          //rest mass temperature of electron
	 nuc=ee*kbperp/me/cc/2./PI,                                    //cyclotron frequency
	 nus=2./9.*nuc*The*The,                                        //peak emissivity frequency
	 Xe=nu/nus;                                                    //ratio of observed frequency to peak emissivity frequency


//RG:nth emissivities -> try to move upwards, figure out dependencies
// jIc_nth = emissivity_coeff * kbperp * jIc_nth_lookuptab;
// jQc_nth = emissivity_coeff * kbperp * jQc_nth_lookuptab;
// jVc_nth = fljVc * (ee*ee*ee*rho*rgrav/2.)/(sqrt(3.)*PI*cc*cc*me) * kbpar * jVc_nth_lookuptab;//plus sign for the "k" vector defined above
// //jVc_nth = jVc; // RG:FIXME SETTING V emissivity equal to old value

//printf("[evalpointzero.cpp]: Use approximate emissivities only to get absorptivities via Kirchhoff's law\n");

jIc_approx=rho*rgrav/2.*sqrt(2.)*PI*ee*ee*nus/6./The/The/cc*Xe*exp(-pow(Xe,(doub)0.3333));       //simple approximation for total emissivity
doub tem=(sqrt(Xe)+pow(2.,11./12.)*pow(Xe,(doub)0.166666));                                //auxiliary
jIc_approx=rho*rgrav/2.*sqrt(2.)*PI*ee*ee*nus/6./The/The/cc*tem*tem*exp(-pow(Xe,(doub)0.3333));  //more accurate approximation for total emissivity from Leung et al. (2011)
aIc_approx=jIc_approx/Bnu;                                                                             //correspondent absorptivity (this is used/stored in ff[4])



// TESTING NTH
if (nth) {
//if (false) {
  //jIc = 0.;                                 // zero out total I   emissivities for testing purposes
  //jIc_nth = 0.;                                // zero out total I   emissivities for testing purposes
  //aIc = 0.; aIc_nth = 0.;                       // zero out total I   absorptivies for testing purposes
  //jQc = 0.; jVc = 0.; jQc_nth = 0.; jVc_nth = 0.; // zero out polarized emissivities for testing purposes
  //aQc = 0.; aVc = 0.; aQc_nth = 0.; aVc_nth = 0.; // zero out polarized absorptivities for testing purposes

  //  printf("jI_old=%e jI_new=%e ratio=%e jI_approximate=%e\n",jIc,jIc_nth,jIc_nth/jIc,jIc_approx);

}




if(jIc!=jIc || jQc!=jQc || jVc!=jVc ){ //check for NANs
	printf("Uncaught otherwise error in radiative transfer...\n Exiting");
	exit(-1);
};
