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
	 sinth,            //sin(polar angle)
	 cossq,            //cos^2(polar angle)
	 sinsq,            //sin^2(polar angle) - to make the code nicer
	 rsq,              //rr*rr
	 Del,              //rsq-2.*rr+asq
	 rhosq;            //rsq+asq*cossq
doub rest[11],         //single record in a fluid simulation file
	 rho,              //physical number density
	 r_T_u[rlen][dd];  //for reading internal energy density radial profile from file
doub iKS[4][4],        //-,+,+,+ signature covariant (upper indices) Kerr-Schild metric
     iBL[4][4],        //-,+,+,+ signature covariant (upper indices) Kerr metric in Boyer-Lindquist coordinates
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
cout << YELLOW"[init.cpp]:"RESET" location of fluid simulation dumps"+fdir+"\n";

//RG:    vvvvvvv HARDWIRE-WARNING, should be a parameter
mintim=1-0.00005*dtimdf*fdiff;       //minimum and maximum proper time at infinity, when the simulation is still evolved together with light ray propagation
maxtim=1+0.00005*dtimdf*fdiff;

//reading initialization files, done once
if(!inited){
    //single record of fluid simulation dump files consists of: rho, u, -u^t, -T^r_t/(rho u^r), u^t, v^r, v^theta, v^phi, B^r, B^theta, B^phi 
	
  cout << YELLOW"[init.cpp]:"RESET" dir="+dir+",astr[sp]="+astr[sp]+",xstr="+xstr+" :"+dir+astr[sp]+xstr+"Tsmap"+astr[sp]+".dat\n";
	//reading mean density & temperature radial profiles
    //RG:
    cout << YELLOW"[init.cpp]:"RESET"location of Tsmap file: "+dir+xstr+"Tsmap"+astr[sp]+".dat\n";

	ifstream r_T_u_file ((dir+astr[sp]+xstr+"Tsmap"+astr[sp]+".dat").c_str(), ios::in);
	for(k=0;k<rlen;k++)
		for(i=0;i<dd;i++)
          r_T_u_file>>r_T_u[k][i];
	r_T_u_file.close();

	//allocating memory for fluid simulation snapshots - takes some time
	for(j=0;j<2*fdiff+1;j++)
		loaded[j]=0;

    // RG: We have done that already in ASTRORAY_main.cpp
    // RG: As long as fdiff is defined globally not locally within cases we can remove the below two-liner
	// for(j=0;j<2*fdiff+1;j++)
	//  	uu[j]=(uuarr) new float[phlen][thlen][rlen][wdd];
    // RG: new should have a matching delete or we leak memory

	//reading coordinate matrix and coordinate transformation matrix
    cout<<YELLOW"[init.cpp]:"RESET" READING dxdxp.dat FILE FROM "<<dir+astr[sp]+xstr+"dxdxp.dat"<<endl;
	ifstream dxp((dir+astr[sp]+xstr+"dxdxp.dat").c_str(), ios::in|ios::binary);
	pbuf=dxp.rdbuf();
	dxp.read(reinterpret_cast<char *>(coord), ndd*thlen*2*sizeof(float));
	dxp.read(reinterpret_cast<char *>(dxdxp), ndd*thlen*4*4*sizeof(float));

    //RG:FIXME TESTING THICKDISK7 dxdxp.dat
    //ndd=650; thlen=128;
	// ifstream dxp((dir+astr[sp]+xstr+"dxdxp-thickdisk7.dat").c_str(), ios::in|ios::binary);
	// pbuf=dxp.rdbuf();
	// dxp.read(reinterpret_cast<char *>(coord), 650*128*2*sizeof(float));
	// dxp.read(reinterpret_cast<char *>(dxdxp), 650*128*4*4*sizeof(float));

	dxp.close();
	
	//polar angles of a "deformed" grid
	for(k=0;k<ndd;k++)
	//for(k=0;k<650;k++)

      //RG: REPORT COORDINATES
      //printf(YELLOW"[init.cpp]:"RESET" r[%d]=%e\n",k,exp(coord[k][0][0]));

		for(i=0;i<thlen;i++) {

          //RG: CHECK FOR NAN
          for(m=0; m<4; m++) {
            if (dxdxp[k][i][1][m]!=dxdxp[k][i][1][m]) {
                cout << YELLOW"[init.cpp]:"RESET" dxdxp coordinates contain nan. Unwise to continue..." << endl;
                exit(1);}
            }
          if (coord[k][i][1]!=coord[k][i][1]) {
            cout << YELLOW"[init.cpp]:"RESET" Theta coordinates contain nan. Unwise to continue..." << endl;
            exit(1);}

			theta[k][i]=coord[k][i][1];
		};

	//radial grid & mean internal energy densities = "temperatures"
	for(k=0;k<ncut;k++) {
        if (coord[k][0][0]!=coord[k][0][0]) { // CHECK FOR NAN
          cout << YELLOW"[init.cpp]:"RESET" Radial coordinates contain nan. Unwise to continue..." << endl;
          exit(1);}

		rtab[k]=exp(coord[k][0][0]);
		Tstab[k]=r_T_u[k][1];

        //RG: REPORT COORDINATES
        //printf(YELLOW"[init.cpp]:"RESET" r[%d]=%e\n",k,rtab[k]);

	};
	maxT=Tstab[0];//maximum internal energy density
	rmin=r_T_u[0][0];//minimum radius (inside the event horizon)

    //reading emissivities for [THERMAL] jI, jQ, jV as well as Faraday effects rQ (Faraday conversion coefficient) and rV (Faraday rotation coefficient)
    cout << YELLOW"[init.cpp]:"RESET" READING lookup???.dat files from "+ASTRORAY_PATH+"/analysis/\n";

	ifstream fI ((ASTRORAY_PATH+"/analysis/lookupjIy.dat").c_str(), ios::in);
	ifstream fQ ((ASTRORAY_PATH+"/analysis/lookupjQy.dat").c_str(), ios::in);
	ifstream fV ((ASTRORAY_PATH+"/analysis/lookupjVy.dat").c_str(), ios::in);
	ifstream frQ ((ASTRORAY_PATH+"/analysis/lookuprQa.dat").c_str(), ios::in);
	ifstream frV ((ASTRORAY_PATH+"/analysis/lookuprVa.dat").c_str(), ios::in);

	for(k=0;k<=Tlen;k++) 
		for(j=0;j<=nWlen;j++) {
			fI>>jI[k][j];
			fQ>>jQ[k][j];
			fV>>jV[k][j];

            //doub nu_Omega_nth = 12000.*pow(1.25,-nWlen_nth/2+j);
            doub nu_Omega_th = 12000.*pow(1.1,-nWlen/2+j);
            doub Theta_e = 0.4*pow(1.05,k);

            // if (nu_Omega_th==12000. && k==0) {
            //   printf(YELLOW"[init.cpp]:"RED" %e %e exp(jI_old)= %e j=%d k=%d #thtable"RESET"\n",nu_Omega_th,Theta_e,exp(jI[k][j]),j,k);
            //   printf(YELLOW"[init.cpp]:"RED" %e %e exp(jQ_old)= %e j=%d k=%d #thtable"RESET"\n",nu_Omega_th,Theta_e,exp(jQ[k][j]),j,k);
            //   printf(YELLOW"[init.cpp]:"RED" %e %e exp(jV_old)= %e j=%d k=%d #thtable"RESET"\n",nu_Omega_th,Theta_e,exp(jV[k][j]),j,k);
            //   //printf(YELLOW"[init.cpp]:"RED" %e %e jI_old=%e j=%d k=%d #thtable"RESET"\n",nu_Omega_th,Theta_e,jI[k][j],j,k);
            // }

		} // for(j=0;j<=nWlen;j++) {

	for(k=0;k<=2*Tlen;k++) {
		frQ>>rQ[k];
		frV>>rV[k];
	}
	fI.close();fQ.close();fV.close();frQ.close();frV.close();


    if (nth) {
      //  stringstream dir_nth="/home/rgold/rt/newthermaltabs/";
      const string       dir_nth="/home/rgold/rt/newthermaltabs/";
      // const string       dir_nth="/home/rgold/codes/rt-git/newthermaltabs/";
      //char dir_nth[64]="/home/rgold/rt/newthermaltabs/";
      const string jI_nth_file=dir_nth+"origlookup_thermalIjcSs_log.dat";

      cout << YELLOW"[init.cpp]:"RESET" READING NEW THERMAL lookup ["+jI_nth_file+",...] but store them in NON-THERMAL origlookup_thermaljcSs_log.dat files from dir="+dir_nth+"\n";

      //RG:TODO CHECK WHETHER FILE EXISTS

    // thermal tabs v1 (sanity check: expect agreement with old tables)
    ifstream fI_nth ((jI_nth_file).c_str(), ios::in);
    //ifstream fI_nth ((dir_nth+"origlookup_thermalIjcSs_log.dat").c_str(), ios::in);
	ifstream fQ_nth ((dir_nth+"origlookup_thermalQjcSs_log.dat").c_str(), ios::in);
	ifstream fV_nth ((dir_nth+"origlookup_thermalVjcSs_log.dat").c_str(), ios::in);

	ifstream faI_nth ((dir_nth+"origlookup_thermalaIT_log.dat").c_str(), ios::in);
	ifstream faQ_nth ((dir_nth+"origlookup_thermalaQT_log.dat").c_str(), ios::in);
	ifstream faV_nth ((dir_nth+"origlookup_thermalaVT_log.dat").c_str(), ios::in);

    //reading emissivities for [NON-THERMAL] jI, jQ, jV as well as Faraday effects rQ (Faraday conversion coefficient) and rV (Faraday rotation coefficient)
    //char nth_power[8]="2.4";
    // cout << "READING in NON-THERMAL origlookup_"<<nth_power<<"jcSs_log.dat files from dir="+dir+"\n";

	//ifstream fI_nth ((dir+"origlookup_"+nth_power+"_IjcSs_log.dat").c_str(), ios::in);
	// ifstream fQ_nth ((dir+"origlookup_"+nth_power+"_QjcSs_log.dat").c_str(), ios::in);
	// ifstream fV_nth ((dir+"origlookup_"+nth_power+"_VjcSs_log.dat").c_str(), ios::in);
    
    // RG: Te rotativities are just fake, we zero them out afterwards
    // RG: problems when lookup sizes differ between th vs nth?
	// ifstream frQ_nth ((dir+"origlookup_thermalrQT_log.dat").c_str(), ios::in);
	ifstream frV_nth ((dir+"origlookup_thermalrVT_log.dat").c_str(), ios::in); 
	ifstream frQ_nth ((dir+"origlookup_thermalrQT_lin.dat").c_str(), ios::in);
	//ifstream frV_nth ((dir+"origlookup_thermalrVT_lin.dat").c_str(), ios::in); 
    // for(int col=0;col<3;col++) 
    // int col=0;
    
      for(k=0;k<=Tlen_nth;k++) 
        for(j=0;j<=nWlen_nth;j++) 

          {
            doub nu_Omega_nth = 12000.*pow(1.25,-nWlen_nth/2+j);
            // doub nu_Omega_th = 12000.*pow(1.1,-nWlen/2+j);
            doub Theta_e = 0.4*pow(1.05,k);

            // READ-IN LOOKUP TABLE FILES
            fI_nth>>jI_nth[k][j];
            fQ_nth>>jQ_nth[k][j];
            fV_nth>>jV_nth[k][j];

            // TRANSLATE from JON's NEW TABLES to ROMAN S. OLD TABLES
            doub conversion=4.44e-30/nu_Omega_nth;
            jI_nth[k][j]=log(exp(jI_nth[k][j])/conversion);
            jQ_nth[k][j]=log(exp(jQ_nth[k][j])/conversion);
            jV_nth[k][j]=log(exp(jV_nth[k][j])/conversion *4./3.);
            
            // if (nu_Omega_nth==12000. && k==0) {// j=60,k=0
            //   printf(YELLOW"[init.cpp]:"RED" %e %e exp(jI_new)= %e j=%d k=%d #nthtable "RESET"\n",nu_Omega_nth,Theta_e,exp(jI_nth[k][j]),j,k);
            //   printf(YELLOW"[init.cpp]:"RED" %e %e exp(jQ_new)= %e j=%d k=%d #nthtable "RESET"\n",nu_Omega_nth,Theta_e,exp(jQ_nth[k][j]),j,k);
            //   printf(YELLOW"[init.cpp]:"RED" %e %e exp(jV_new)= %e j=%d k=%d #nthtable "RESET"\n",nu_Omega_nth,Theta_e,exp(jV_nth[k][j]),j,k);
            // }

            faI_nth>>aI_nth[k][j];
            faQ_nth>>aQ_nth[k][j];
            faV_nth>>aV_nth[k][j];
            aI_nth[k][j]=log(exp(aI_nth[k][j])/conversion);
            aQ_nth[k][j]=log(exp(aQ_nth[k][j])/conversion);
            aV_nth[k][j]=log(exp(aV_nth[k][j])/conversion*4./3.);
		}

	for(k=0;k<=2*Tlen_nth;k++) {

			frQ_nth>>rQ_nth[k];
			frV_nth>>rV_nth[k];

	}
	fI_nth.close();fQ_nth.close();fV_nth.close();frQ_nth.close();frV_nth.close();
    



	inited=true; //repeat once //RG:never repeated? if !(inited){}
    };
    
    if (false){ // CHECK lookup tables
    int inspect=0;
    int inspect_W=nWlen/2,inspect_T=0;
    int inspect_W_nth=nWlen_nth/2,inspect_T_nth=0;
    printf(YELLOW"[init.cpp]:"RESET" frequency=%E temperature=%E\n",12000.*pow(1.1, -nWlen/2+(inspect_W/2)),0.4*pow(1.05, inspect_T));
    printf(YELLOW"[init.cpp]:"RESET" frequency_nth=%E temperature_nth=%E\n",12000.*pow(1.25, -nWlen_nth/2+(inspect_W_nth/2)),0.4*pow(1.05, inspect_T_nth));
    printf(YELLOW"[init.cpp]:"RESET" exp(jI[%d][%d])=%E exp(jI_nth[%d][%d])=%E\n",inspect_T,inspect_W,exp(jI[inspect][inspect]),inspect_T_nth,inspect_W_nth,exp(jI_nth[inspect][inspect]));
    }
};

//reading fluid simulation dump files from disk
for(i=-fdiff;i<=fdiff;i++) {
	int ioff=i+fdiff;         //consecutive number in an array of loaded dump files
	int fx=fnum+i;            //ID of a dump file

	//checking if simulation dump/snapshot of interest is already in memory
	bool reuse=false;         
	for(j=0;j<2*fdiff+1;j++) 
		if(loaded[j]==fx) {
			printf(YELLOW"[init.cpp]:"RESET"Found coincidence for fnum=%d\n",fx);
			uuarr ptr=uu[ioff]; uu[ioff]=uu[j]; uu[j]=ptr;
			loaded[ioff]=fx;
			reuse=true;
			break;
		};
	if(reuse)continue;//the snapshot we want to load is already in memory
    
	//reading the simulation snapshot
    //RG:TODO HAVE ENCLOSING if(MODEL==SIMULATION) STATEMENT (AS OPPOSED TO "ANALYTIC")
	stringstream sstr;
	sstr<<setfill('0')<<setw(4)<<fx;
	ifstream fline((fdir+"fieldline"+sstr.str()+".bin").c_str(),ios::in|ios::binary);
	pbuf=fline.rdbuf();
	int fsize=pbuf->pubseekoff (0,ios::end),                //file size
		tosize=wdd*phlen*thlen*rlen*sizeof(float);          //size of binary region
	pbuf->pubseekpos(fsize-tosize);                         //discarding the text region in the beginning of dump file

    /* RG: TESTING remove header in file */
	// pbuf->pubseekpos(0);
    //fline.read(reinterpret_cast<char *>(uu[ioff]), fsize-tosize);// -> 0
    //char header[156];
    //FILE *my_fp=fopen("fieldline1000.bin","r");
    //fscanf(my_fp,156);// -> 0
    //my_fp.close();
    //cout << "156 header bytes: " << fline << endl;
    //exit(1);
    /***************/

	fline.read(reinterpret_cast<char *>(uu[ioff]), tosize); //reading as binary
    if(fline.good())
		printf(YELLOW"[init.cpp]:"RESET" fieldline %d read\n",fx);
	else { 
      printf(YELLOW"[init.cpp]:"RED" Something is wrong loading fluid simulation dumps\n...EXITING...\n"RESET);
		exit(-1);
	};
	fline.close();
	loaded[ioff]=fx;
};



        /********************************************************
         * OVERWRITE DATA WITH AVERY'S FORCE-FREE TOY JET MODEL *
         ********************************************************/

        if ( !strcmp(avery_toy_jet,"yes") ){
          if (fdiff!=0) {
            fdiff=0;
            printf(YELLOW"[init.cpp]:"RESET"You called Avery's toyjet model with fdiff=%d. Will set fdiff=0 (model is stationary).\n",fdiff);
          }
          //RG: SETUP AVERY'S TOYJET + RIAF MODEL
          // consider use of PREPROCESSOR DIRECTIVE: #ifdef TOYJET ?
          #include "setup_averys_toyjet.cpp"
        }



/************************************
 * RG: BCs near the poles ...WIP... *
 ************************************/

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
/**********************************************************************************************************/



rcut=r_T_u[ncut-1][0];                                     //define radius up to which the simulation converged
rhopo=-log(rhoout/rhonor/r_T_u[ncut-1][2])/log(rrmax/rcut);//slope of powerlaw density (extension) profile
Upo=-log(Tout/r_T_u[ncut-1][1])/log(rrmax/rcut);           //slope of powerlaw temperature (extension) profile

//----------------tests------------------------
//testing sensitivity to slopes of power-law extensions of density, temperature, and magnetic field
	if(echeck1)rhopo+=0.2;
	if(echeck2)Upo-=0.1;
	if(echeck3)Bpo+=0.2;
//-------------end-of-tests--------------------

rhocon=rhonor*r_T_u[ncut-1][2];   //physical density at radial convergence boundary
Ucon=r_T_u[ncut-1][1];            //temperature at radial convergence boundary
Bnor=sqrt(4.*PI*rhonor*mp)*cc; //normalization of physical magnetic field, negative normalization is also allowed

//extending radial temperature profile outwards as a powerlaw
for(k=ncut;k<ndd;k++){
	rtab[k]=rcut*pow((doub)rrmax/rcut,(k-ncut+(doub)1.)/(doub)(ndd-ncut));
	Tstab[k]=Ucon*pow((doub)rtab[k]/rcut,(doub)-Upo);
};

int nx=60;                       //radial point where accretion rate is computed
doub rx=exp(coord[nx-1][0][0]),  //radius at nx
	 u0;                         //u^t
printf(YELLOW"[init.cpp]:"RESET" HARDWIRED r[nx=%d]=%e (radius where mdot is evaluated)\n",nx,rx);
rate=0.;

//computing the accretion rate in code units
for(k=0;k<thlen-1;k++)
	for(i=0;i<phlen;i++)
      for(j=0;j<=2*fdiff;j++) { // RG: j -> t ?
			u0=(*uu[j])[i][k][nx-1][4];
			usp[j][i][k][0]=u0;
			usp[j][i][k][1]=u0*(*uu[j])[i][k][nx-1][5];//u^r
			usp[j][i][k][2]=u0*(*uu[j])[i][k][nx-1][6];//u^theta
			usp[j][i][k][3]=u0*(*uu[j])[i][k][nx-1][7];//u^phi
			uspKS[j][i][k]=0.;
			for(m=0;m<4;m++)
				uspKS[j][i][k]+=dxdxp[nx-1][k][1][m]*usp[j][i][k][m];//radial velocity in KS metric
			rate+=2./phlen*PI*(rx*rx+ asq*theta[nx-1][k]*theta[nx-1][k+1])*(*uu[j])[i][k][nx-1][0]*uspKS[j][i][k]*(theta[nx-1][k]-theta[nx-1][k+1]);

            //RG: Check for NAN
            if (isnan(rate)){//!=rate){
              printf(YELLOW"[init.cpp]:"RESET"rate=%e,i=%d,j=%d,k=%d,rhonor=%e,rgrav=%e,fdiff=%d\n",rate,i,j,k,rhonor,rgrav,fdiff);
              printf(YELLOW"[init.cpp]:"RESET"rate=%e,i=%d,j=%d,k=%d,rx=%e,asq=%e,theta[nx-1][%d]=%e\n",rate,i,j,k,rx,asq,k+1,theta[nx-1][k+1]);
              printf(YELLOW"[init.cpp]:"RESET"rate=%e,uu[%d][%d][%d][nx-1][0]=%e,uspKS[j][i][k]=%e\n",rate,j,i,k,(*uu[j])[i][k][nx-1][0],uspKS[j][i][k]);
              printf(YELLOW"[init.cpp]:"RESET"rate=%e,dxdxp[nx-1][k][1][0]=%e,usp[j][i][k][0]=%e\n",rate,dxdxp[nx-1][k][1][0],usp[j][i][k][0]); 
              printf(YELLOW"[init.cpp]:"RESET"rate=%e,dxdxp[nx-1][k][1][1]=%e,usp[j][i][k][1]=%e\n",rate,dxdxp[nx-1][k][1][1],usp[j][i][k][1]);
              printf(YELLOW"[init.cpp]:"RESET"rate=%e,dxdxp[nx-1][k][1][2]=%e,usp[j][i][k][2]=%e\n",rate,dxdxp[nx-1][k][1][2],usp[j][i][k][2]); 
              printf(YELLOW"[init.cpp]:"RESET"rate=%e,dxdxp[nx-1][k][1][3]=%e,usp[j][i][k][3]=%e\n",rate,dxdxp[nx-1][k][1][3],usp[j][i][k][3]); 

              exit(1); /***************************************** EXIT *****************/
            }

		};

rate*=rhonor*rgrav*rgrav*cc*mp/(2*fdiff+1);//accretion rate in physical units

//initializing Te(Ts)+Tp(Ts) solver, standard for GSL
cout << YELLOW"[init.cpp]:"RESET" Electron Temperature Solver..." << endl;

doub acc=2e-3,   //relative accuracy //RG: SHOULD BE (CONSTANT) GLOBAL, USER SHOULD NOT HAVE TO HUNT THIS DOWN IN THE CODE
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
step=-0.001*rz; //initial step //RG:HARDWIRED
ts[0]=IT[0];tp[0]=IT[0];te[0]=IT[0];

//actual temperature solver
while (rz > 1.001*rmin) {    //while outside of minimum radius (inside of event horizon). 
	stNt++; 
	if(-step>0.008*rz)       
		step=-0.005*rz;      //artificially limiting the step
	int status = gsl_odeiv_evolve_apply (ez, cz, sz, &sysT, &rz, 1.001*rmin, &step, IT);
	if(status!=0){
		printf(YELLOW"[init.cpp]:"RED"Temperature solver error\n Exiting"RESET);
		exit(-1);
	};
	te[stNt]=IT[0];
	tp[stNt]=IT[1];
	if((rz<6) && (fl)) {    //catching electron temperature & Tp/Te ratio at 6M distance from the BH. Step size is so small that we don't iterate for 6M precisely
		fl=false;
		TpTe=tp[stNt]/te[stNt];
		Te6=te[stNt];
		printf(YELLOW"[init.cpp]:"RESET" fn=%d stN=%d r=%.2fM Tp/Te=%.2f Te=%.3e rate=%.3e he=%.3f rho=%.3e\n", fnum, stNt,rz, TpTe, Te6, rate*year/Msun,heat,rhonor);
	};
;};
maxT=ts[stNt];minT=ts[0];

if(stNt>maxst){            //check if solver exceeded the number of steps
	printf(YELLOW"[init.cpp]:"RED"Temperature solver error. Exceeded maximum steps allowed! \nEXITING\n"RESET);
	exit(-1);
};
gsl_odeiv_evolve_free (ez);//free memory
gsl_odeiv_control_free (cz);
gsl_odeiv_step_free (sz);



/*************************
 * OUTER RADIAL BOUNDARY *
 *************************/

//compute quantities in the co-moving frame at the convergence boundary (given by rr=rcut=const)
xrn=ncut-1; //point just inside convergence boundary
rr=rcut;    

/****************************************************************************/
for(nt=0;nt<thlen;nt++){ //for all POLAR angles at the convergence boundary *
	costh=theta[xrn][nt];
	cossq=costh*costh;
	sinsq=1-cossq;
    sinth=sqrt(sinsq);
	rsq=rr*rr;
	rhosq=rsq+asq*cossq;
	Del=rsq-2.*rr+asq; // Lambda_squared in http://grwiki.physics.ncsu.edu/wiki/Kerr_Black_Hole
	doub temp=2.*rr/rhosq;

	// contravariant BL metric matrix
	// iBL[0][0]=-(Del/rhosq) - a*sinsq/rhosq;
	// iBL[0][1]=0.;
    // iBL[0][2]=0.;
	// iBL[0][3]=Del/rhosq * a *sinsq - sinsq/rhosq*(rsq+asq)*a;
	// iBL[1][0]=iBL[0][1];
	// iBL[1][1]=rhosq/Del;
	// iBL[1][2]=0.;iBL[1][3]=0.;
	// iBL[2][0]=iBL[0][2];
	// iBL[2][1]=iBL[1][2];
	// iBL[2][2]=rhosq;
	// iBL[2][3]=0.;
	// iBL[3][0]=iBL[0][3];
	// iBL[3][1]=iBL[1][3];
	// iBL[3][2]=iBL[2][3];
	// iBL[3][3]=Del/rhosq*a*sinsq + sinsq/rhosq*(rsq+asq);

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
	
    /****************************************************/
	for(np=0;np<phlen;np++){ //for all AZIMUTHAL angles *

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


    }; 	// for(np=0;np<phlen;np++){ //for all AZIMUTHAL angles *
 }; // for(nt=0;nt<thlen;nt++){ //for all POLAR angles at the convergence boundary *



printf("\n"YELLOW"[init.cpp]:"RESET" Exiting [init.cpp]\n");

return 0; // init()
}


