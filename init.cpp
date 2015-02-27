int init(int sp, int fmin, int fmax, int sep) { //RG: sep is never used here... remove!
	int m, k, j, n, i,
		xrn,           //radial index at the boundary of convergence
		nt,            //theta index
		np;            //phi index
bool fl=true,          //for catching 6M distance from BH in temperature calculation
	 inited=false;     //whether initialization already performed
filebuf *pbuf;         //buffer for reading file with offset

doub rr,               //radius
	 costh,            //cos(polar angle)
	 cossq,            //cos^2(polar angle)
	 sinsq,            //sin^2(polar angle) - to make the code nicer
	 rsq,              //rr*rr
	 Del,              //rsq-2.*rr+asq
	 rhosq;            //rsq+asq*cossq
doub rest[11],         //single record in a fluid simulation file
	 rho,              //physical number density
	 xx[rlen][dd];     //for reading internal energy density radial profile from file
doub iKS[4][4],        //-,+,+,+ signature covariant (upper indices) Kerr-Schild metric
	 iMKS[4][4],       //modified Kerr-Schild metric (code coordinates)
	 MKStoKS[4][4],    //(covariant) vector transformation from MKS to KS
	 u[4],             //4-velocity in MKS
	 ulo[4],           //covariant (low indices) 4-velocity in MKS
	 uKS[4],           //4-velocity in KS
	 uloKS[4],         //covariant (low indices) 4-velocity in KS
	 Bi[4],            //3-vector of magnetic field (Bi[0]==0) (NOT physical magnetic field)
	 Bup[4],           //covariant 4-vector of magnetic field in MKS
	 BupKS[4],         //covariant 4-vector of magnetic field in KS
	 sc[4];            //auxiliary for vector normalization

a=atab[sp]; asq=a*a;                 //define spin value for chosen fluid simulation
ncut=ncuttab[sp];                    //define radial index at the point, where the fluid simulation is barely converged
string fdir=adir+astr[sp]+fieldstr;  //location of fluid simulation dumps
cout << "location of fluid simulation dumps"+fdir+"\n";

//RG:    vvvvvvv HARDWIRE-WARNING, should be a parameter
mintim=1-0.00005*dtimdf*fdiff;       //minimum and maximum proper time at infinity, when the simulation is still evolved together with light ray propagation
maxtim=1+0.00005*dtimdf*fdiff;

//reading initialization files, done once
if(!inited){
    //single record of fluid simulation dump files consists of: rho, u, -u^t, -T^r_t/(rho u^r), u^t, v^r, v^theta, v^phi, B^r, B^theta, B^phi 
	
  cout << "dir="+dir+"\n,astr[sp]="+astr[sp]+"\n,xstr="+xstr+" :"+dir+astr[sp]+xstr+"Tsmap"+astr[sp]+".dat\n";
	//reading mean density & temperature radial profiles
    //RG:
    cout << "location of Tsmap file: "+dir+xstr+"Tsmap"+astr[sp]+".dat\n";

	ifstream faire ((dir+astr[sp]+xstr+"Tsmap"+astr[sp]+".dat").c_str(), ios::in);
	for(k=0;k<rlen;k++)
		for(i=0;i<dd;i++)
			faire>>xx[k][i];
	faire.close();

	//allocating memory for fluid simulation snapshots - takes some time
	for(j=0;j<2*fdiff+1;j++)
		loaded[j]=0;
	for(j=0;j<2*fdiff+1;j++)
		uu[j]=(uuarr) new float[phlen][thlen][rlen][wdd];

	//reading coordinate matrix and coordinate transformation matrix
	ifstream dxp((dir+astr[sp]+xstr+"dxdxp.dat").c_str(), ios::in|ios::binary);
	pbuf=dxp.rdbuf();
	dxp.read(reinterpret_cast<char *>(coord), ndd*thlen*2*sizeof(float));
	dxp.read(reinterpret_cast<char *>(dxdxp), ndd*thlen*4*4*sizeof(float));
	dxp.close();
	
	//polar angles of a "deformed" grid
	for(k=0;k<ndd;k++)
		for(i=0;i<thlen;i++) {
			theta[k][i]=coord[k][i][1];
		};
	//radial grid & mean internal energy densities = "temperatures"
	for(k=0;k<ncut;k++) {
		rtab[k]=exp(coord[k][0][0]);
		Tstab[k]=xx[k][1];
	};
	maxT=Tstab[0];//maximum internal energy density
	rmin=xx[0][0];//minimum radius (inside the event horizon)

    //reading emissivities jI, jQ, jV as well as Faraday effects rQ (Faraday conversion coefficient) and rV (Faraday rotation coefficient)
    cout << "READING in lookup???.dat files from dir="+dir+"\n";

	ifstream fI ((dir+"lookupjIy.dat").c_str(), ios::in);
	ifstream fQ ((dir+"lookupjQy.dat").c_str(), ios::in);
	ifstream fV ((dir+"lookupjVy.dat").c_str(), ios::in);
	ifstream frQ ((dir+"lookuprQa.dat").c_str(), ios::in);
	ifstream frV ((dir+"lookuprVa.dat").c_str(), ios::in);
	for(k=0;k<=Tlen;k++) 
		for(j=0;j<=nWlen;j++) {
			fI>>jI[k][j];
			fQ>>jQ[k][j];
			fV>>jV[k][j];
		}
	for(k=0;k<=2*Tlen;k++) {
		frQ>>rQ[k];
		frV>>rV[k];
	}
	fI.close();fQ.close();fV.close();frQ.close();frV.close();

	inited=true; //repeat once
};

//reading fluid simulation dump files from disk
for(i=-fdiff;i<=fdiff;i++) {
	int ioff=i+fdiff;         //consecutive number in an array of loaded dump files
	int fx=fnum+i;            //ID of a dump file

	//checking whether the simulation dump/snapshot of interest is already in memory
	bool reuse=false;         
	for(j=0;j<2*fdiff+1;j++) 
		if(loaded[j]==fx) {
			printf("Found coincidence for fnum=%d\n",fx);
			uuarr ptr=uu[ioff]; uu[ioff]=uu[j]; uu[j]=ptr;
			loaded[ioff]=fx;
			reuse=true;
			break;
		};
	if(reuse)continue;//the snapshot we want to load is already in memory
    
	//reading the simulation snapshot
	stringstream sstr;
	sstr<<setfill('0')<<setw(4)<<fx;
	ifstream fline((fdir+"fieldline"+sstr.str()+".bin").c_str(),ios::in|ios::binary);
	pbuf=fline.rdbuf();
	int fsize=pbuf->pubseekoff (0,ios::end),                //file size
		tosize=wdd*phlen*thlen*rlen*sizeof(float);          //size of binary region
	pbuf->pubseekpos(fsize-tosize);                         //discarding the text region in the beginning of dump file
	fline.read(reinterpret_cast<char *>(uu[ioff]), tosize); //reading as binary
    if(fline.good())
		printf("fieldline %d read\n",fx);
	else { 
		printf("Something is wrong loading fluid simulation dumps \n Exiting");
		exit(-1);
	};
	fline.close();
	loaded[ioff]=fx;
};

/**********************************************************************************************************/
//RG: MODIFY FLUID DATA HERE?
/**********************************************************************************************************/
// single fluid simulation dump files contains: 
// uu[?]:     [0    1   2     3                4    5    6        7      8    9        10   ]
//             rho, u, -u^t, -T^r_t/(rho u^r), u^t, v^r, v^theta, v^phi, B^r, B^theta, B^phi
//for(theta_index=0;theta_index<thlen-1;theta_index++)
// int theta_bc_cells = 5;
// //RG: HMM...                        vv
//  for(theta_index=0;theta_index<thlen-1;theta_index++) // only boundaries
//    for(phi_index=0;phi_index<phlen;phi_index++)
// 	  for(r_index=0;r_index<rlen;r_index++)
// 		for(snapshot_iteration_index=0;snapshot_iteration_index<=2*fdiff;snapshot_iteration_index=0++) {

//           if (theta_index <= theta_bc_cells){
//             // extrapolate rho,ug,v1,v3,b1,b3 as symmetric (extrapolate as just constant)
//             // v2,b2 as anti-symmetric (passed through zero at pole, so linear from cell 6 to pole where has value of zero)
//           }
//           else if (theta_index > thlen-theta_bc_cells){
//             // extrapolate rho,ug,v1,v3,b1,b3 as symmetric (extrapolate as just constant)
//             // v2,b2 as anti-symmetric (passed through zero at pole, so linear from cell 6 to pole where has value of zero)
//           }
//           //RG: FIXME: careful about T units (could be in K)
//           // doub Y_e=0.62; doub mubar=0.62; // solar abundances
//           // doub T_e=6.; // =3e10K (T_rest_electron=1=me c^2,me=511kev)
//           // doub rho_orig = uu[snapshot_iteration_index][phi_index][theta_index][r_index][0]; // grmhd (HARM) code
//           // doub u_orig = uu[snapshot_iteration_index][phi_index][theta_index][r_index][1]; // grmhd (HARM) code
//           // doub gamma_orig=4./3.; // 13./9.; // u_orig/P+1.;
//           // //RG:                               vv
//           // doub u_Te_cap = Y_e*rho_orig/(mubar*mb)*kb*T_e/(gamma_orig-1.); // = n_e kb T_e/(gamma-1) 
//           // // set u_g=u_Te_cap get T_e the desired isothermal (constant) value, e.g. T_e = 3\times 10^{10}K.;
//           // //RG:FIXME magn_cap not available at this stage in the code... use cmdline argument directly?
//           // //RG:FIXME      vvv
//           // // doub magn     = ???;// see [evalpointzero.cpp]
//           // // doub magn_cap = 4; // see [evalpointzero.cpp]
//           // //(*uu[snapshot_iteration_index])[i][theta_index][r_index][0] = u_orig*(1.-exp(-magn/magn_cap)) + u_Te_cap*(exp(-magn/magn_cap));
//           (*uu[snapshot_iteration_index])[i][theta_index][r_index][0] = u_orig*(1.-exp(-magn/magn_cap)) + u_Te_cap*(exp(-magn/magn_cap));
//         }
/**********************************************************************************************************/
/**********************************************************************************************************/

rcut=xx[ncut-1][0];                                     //define radius up to which the simulation converged
rhopo=-log(rhoout/rhonor/xx[ncut-1][2])/log(rrmax/rcut);//slope of powerlaw density (extension) profile
Upo=-log(Tout/xx[ncut-1][1])/log(rrmax/rcut);           //slope of powerlaw temperature (extension) profile

//----------------tests------------------------
//testing sensitivity to slopes of power-law extensions of density, temperature, and magnetic field
	if(echeck1)rhopo+=0.2;
	if(echeck2)Upo-=0.1;
	if(echeck3)Bpo+=0.2;
//-------------end-of-tests--------------------

rhocon=rhonor*xx[ncut-1][2];   //physical density at radial convergence boundary
Ucon=xx[ncut-1][1];            //temperature at radial convergence boundary
Bnor=sqrt(4.*PI*rhonor*mp)*cc; //normalization of physical magnetic field, negative normalization is also allowed

//extending radial temperature profile outwards as a powerlaw
for(k=ncut;k<ndd;k++){
	rtab[k]=rcut*pow((doub)rrmax/rcut,(k-ncut+(doub)1.)/(doub)(ndd-ncut));
	Tstab[k]=Ucon*pow((doub)rtab[k]/rcut,(doub)-Upo);
};

int nx=60;                       //radial point where accretion rate is computed
doub rx=exp(coord[nx-1][0][0]),  //radius at nx
	 u0;                         //u^t
rate=0.;

//computing the accretion rate in code units
for(k=0;k<thlen-1;k++)
	for(i=0;i<phlen;i++)
		for(j=0;j<=2*fdiff;j++) {
			u0=(*uu[j])[i][k][nx-1][4];
			usp[j][i][k][0]=u0;
			usp[j][i][k][1]=u0*(*uu[j])[i][k][nx-1][5];//u^r
			usp[j][i][k][2]=u0*(*uu[j])[i][k][nx-1][6];//u^theta
			usp[j][i][k][3]=u0*(*uu[j])[i][k][nx-1][7];//u^phi
			uspKS[j][i][k]=0.;
			for(m=0;m<4;m++)
				uspKS[j][i][k]+=dxdxp[nx-1][k][1][m]*usp[j][i][k][m];//radial velocity in KS metric
			rate+=2./phlen*PI*(rx*rx+ asq*theta[nx-1][k]*theta[nx-1][k+1])*(*uu[j])[i][k][nx-1][0]*uspKS[j][i][k]*(theta[nx-1][k]-theta[nx-1][k+1]);
		};
rate*=rhonor*rgrav*rgrav*cc*mp/(2*fdiff+1);//accretion rate in physical units

//initializing Te(Ts)+Tp(Ts) solver, standard for GSL
doub acc=2e-3,   //relative accuracy
	 IT[2],      //ODE vector (Te+Tp)
	 step,       //step
	 rz;         //radius
gsl_odeiv_step *sz; 
gsl_odeiv_control *cz; 
gsl_odeiv_evolve *ez;
const gsl_odeiv_step_type *Tz; 
Tz=gsl_odeiv_step_rk2;//Runge-Kutta 2-nd order
gsl_odeiv_system sysT = {solvetemperature, NULL, 2,NULL};
ez = gsl_odeiv_evolve_alloc(2);
cz = gsl_odeiv_control_standard_new(0.0, acc, 1.0, 0.0);
sz = gsl_odeiv_step_alloc (Tz, 2);

//setting up the first step of integration inwards from the outer boundary
IT[0]=Tstab[ndd-1];
IT[1]=IT[0];
stNt=0;
rz=rrmax;
step=-0.001*rz; //initial step
ts[0]=IT[0];tp[0]=IT[0];te[0]=IT[0];

//actual temperature solver
while (rz > 1.001*rmin) {    //while outside of minimum radius (inside of event horizon). 
	stNt++; 
	if(-step>0.008*rz)       
		step=-0.005*rz;      //artificially limiting the step
	int status = gsl_odeiv_evolve_apply (ez, cz, sz, &sysT, &rz, 1.001*rmin, &step, IT);
	if(status!=0){
		printf("Temperature solver error\n Exiting");
		exit(-1);
	};
	te[stNt]=IT[0];
	tp[stNt]=IT[1];
	if((rz<6) && (fl)) {    //catching electron temperature & Tp/Te ratio at 6M distance from the BH. Step size is so small that we don't iterate for 6M precisely
		fl=false;
		TpTe=tp[stNt]/te[stNt];
		Te6=te[stNt];
		printf("fn=%d stN=%d r=%.2fM Tp/Te=%.2f Te=%.3e rate=%.3e he=%.3f rho=%.3e\n", fnum, stNt,rz, TpTe, Te6, rate*year/Msun,heat,rhonor);
	};
;};
maxT=ts[stNt];minT=ts[0];

if(stNt>maxst){            //check if solver exceeded the number of steps
	printf("Temperature solver error. Exceeded buffer! \n Exiting");
	exit(-1);
};
gsl_odeiv_evolve_free (ez);//free memory
gsl_odeiv_control_free (cz);
gsl_odeiv_step_free (sz);


//compute quantities in the co-moving frame at the convergence boundary (given by rr=rcut=const)
xrn=ncut-1; //point just inside convergence boundary
rr=rcut;    

for(nt=0;nt<thlen;nt++){ //for all polar angles at the convergence boundary
	costh=theta[xrn][nt];
	cossq=costh*costh;
	sinsq=1-cossq;
	rsq=rr*rr;
	rhosq=rsq+asq*cossq;
	Del=rsq-2.*rr+asq;
	doub temp=2.*rr/rhosq;
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

	//define MKS to KS transformation for given coordinates
	for(i=0;i<4;i++)
		for(k=0;k<4;k++)
			MKStoKS[i][k]=dxdxp[xrn][nt][i][k];
	//initialize with zeros metric tensor in MKS
	for(i=0;i<4;i++)
		for(k=0;k<4;k++)
			iMKS[i][k]=0.;
	//compute metric tensor in MKS via a basic tensor transformation
	for(i=0;i<4;i++)
		for(k=0;k<4;k++)
			for(j=0;j<4;j++)
				for(m=0;m<4;m++)
					iMKS[i][k]+=MKStoKS[m][i]*iKS[m][j]*MKStoKS[j][k];
	
	for(np=0;np<phlen;np++){ //for all azimuthal angles
		//read fluid quantities at a given point
		for(m=0;m<11;m++)
			rest[m]=(*uu[fdiff])[np][nt][xrn][m];
		rho=rest[0]*rhonor;              //density
		rest[1]*=mp*cc*cc/3/kb/(rest[0]);//temperature
		u[0]=rest[4];                    //4-velocity
		u[1]=u[0]*rest[5];
		u[2]=u[0]*rest[6];
		u[3]=u[0]*rest[7];
		Bi[1]=Bnor*rest[8];              //3-vector of magnetic field
		Bi[2]=Bnor*rest[9];
		Bi[3]=Bnor*rest[10];
        
		//writing into a collector variable
		uext[np][nt][3]=rho;
        uext[np][nt][4]=rest[1];
		//4-velocity in KS coordinates
		for(i=1;i<4;i++){
			uKS[i]=0.;
			for(k=0;k<4;k++)
				uKS[i]+=MKStoKS[i][k]*u[k];
		};
		//computing uKS[0] precisely. Simple coordinate transformation leads to loss of accuracy
		doub expr=fabs((uKS[1]*iKS[0][1] + uKS[3]*iKS[0][3])* (uKS[1]*iKS[0][1] + uKS[3]*iKS[0][3]) - iKS[0][0]*(1. + uKS[1]*uKS[1]*iKS[1][1] + 2.*uKS[1]*uKS[3]*iKS[1][3] + uKS[2]*uKS[2]*iKS[2][2] + uKS[3]*uKS[3]*iKS[3][3]));
		uKS[0]=-((uKS[1]*iKS[0][1] + uKS[3]*iKS[0][3] + sqrt(expr))/iKS[0][0]);
		//covariant 4-vector of velocity in KS coordinates
		for (m=0;m<4;m++){
			uloKS[m]=0.;
			for(j=0;j<4;j++)
				uloKS[m]+=iKS[m][j]*uKS[j];
		};
		u[0]=uKS[0];//works only in simple MKS, where KS->MKS is a spatial transformation
		//covariant 4-vector of velocity in MKS coordinates
		for(m=0;m<4;m++){
			ulo[m]=0.;
			for(j=0;j<4;j++)
				ulo[m]+=iMKS[m][j]*u[j];
		};
		     //transformation matrix between KS contravariant vector and locally flat co-moving reference frame - requires a correction (see below)
		     //this algorithm doesn't require a correction in BL coordinates. Better algorithm may exist for KS coordinates...
		doub yyKS[4][4]={{uKS[0], uKS[1], uKS[2], uKS[3]},
		                 {uKS[0]*uloKS[1], -(uKS[0]*uloKS[0] + uKS[3]*uloKS[3]), 0, uKS[3]*uloKS[1]},
		                 {uloKS[2]*uKS[0], uloKS[2]*uKS[1],uloKS[2]*uKS[2]+1, uloKS[2]*uKS[3]},
		                 {-uloKS[3]/uloKS[0], 0, 0, 1}}, 
			 //transformation matrix between KS covariant vector and locally flat co-moving reference frame
			 yyloKS[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},
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
		};
		//testing that velocity in locally flat co-moving reference frame is (1,0,0,0)
		doub test[4]={0.,0.,0.,0.};
		for(m=0;m<4;m++){
			for(j=0;j<4;j++)
				test[m]+=yyloKS[m][j]*uKS[j];
		};
		//computing co-moving magnetic field (finally after so much effort...)
		doub B[4]={0.,0.,0.,0.};
		for(m=0;m<4;m++){
			for(j=0;j<4;j++)
				B[m]+=yyloKS[m][j]*BupKS[j];
		};
		//writing into a collector variable
        uext[np][nt][0]=B[1];
        uext[np][nt][1]=B[2];
        uext[np][nt][2]=B[3];
    };
};

printf("\nExiting [init.cpp]\n");

return(0);
}
