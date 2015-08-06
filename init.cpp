int init(int sp, int fmin, int fmax, int sep) { //RG: sep is never used here... remove!
  
    time_t t_b4_init = clock(); // PROFILING

	int m, k, j, n, i,
		r_ext_boundary_idx,           //radial index at the boundary of convergence
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
	 Delta,            //rsq-2.*rr+asq
	 rhosq;            //rsq+asq*cossq //akas Sigma
doub rest[11],         //single record in a fluid simulation file
	 rho,              //physical number density //RG: meaning n (not rho?)
	 r_T_u[rlen][dd];  //for reading internal energy density radial profile from file
doub iKS[4][4],        //-,+,+,+ signature covariant (upper indices) Kerr-Schild metric
     KS[4][4],         //RG: -,+,+,+ signature covariant? (lower indices) Kerr-Schild metric
     BL[4][4],         //-,+,+,+ signature covariant (lower indices) Kerr metric in Boyer-Lindquist coordinates
     iBL[4][4],        //-,+,+,+ signature contravariant/inverse (upper indices) Kerr metric in Boyer-Lindquist coordinates
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

     // AVERY's MODEL
     doub F[2][rlen];      // F(ψ) array -> loading the jet
     doub FofPsi[2][rlen]; // F(ψ) array -> loading the jet (tmp array)
     doub gamma;           // Avery's model
     doub b_F_sq;           // Avery's model


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
	
	
	for(int r_idx=0;r_idx<ndd;r_idx++) // polar angles of a "deformed" grid // RG:?
		for(int th_idx=0;th_idx<thlen;th_idx++) {

          // for(int mu=0; mu<4; mu++) {
          //   if (isnan(dxdxp[r_idx][th_idx][1][mu])) {
          for(int mu=0; mu<4; mu++) for(int nu=0; nu<4; nu++) {

              // RG: Is setting dxdxp=1 inconsistent with rest of code?
              // if ( !strcmp(avery_toy_jet,"yes") ) {
              //   dxdxp[r_idx][th_idx][mu][nu] = 0;
              //   if (mu==mu) dxdxp[r_idx][th_idx][mu][nu] = 1;
              //   if (mu+nu+r_idx+th_idx==0) printf(YELLOW"[init.cpp]: "RED"DEFORMED GRID DEACTIVATED (dxdxp=1)\n"RESET);
              // }

              //RG: CHECK FOR NAN
              if (isnan(dxdxp[r_idx][th_idx][mu][nu])) {
                cout << YELLOW"[init.cpp]:"RESET" dxdxp coordinates contain nan. Unwise to continue..." << endl;
                exit(1);}
            }
          if (isnan(coord[r_idx][th_idx][1])) {
            cout << YELLOW"[init.cpp]:"RESET" Theta coordinates contain nan. Unwise to continue..." << endl;
            exit(1);}

          theta[r_idx][th_idx]=coord[r_idx][th_idx][1]; //RG: theta != costh
		};

	//radial grid & mean internal energy densities = "temperatures"
	for(k=0;k<ncut;k++) {
      if (isnan(coord[k][0][0])) { // CHECK FOR NAN
        cout<<YELLOW"[init.cpp]:"RESET" Radial coordinates contain nan. Unwise to continue..."<<endl;
        exit(1);
      }

      rtab[k]=exp(coord[k][0][0]);
      Tstab[k]=r_T_u[k][1];
	};

	maxT=Tstab[0];//maximum internal energy density
	rmin=r_T_u[0][0];//minimum radius (inside the event horizon)

    /*****************************************************************************/
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
		} // for(j=0;j<=nWlen;j++) {

	for(k=0;k<=2*Tlen;k++) {
		frQ>>rQ[k];
		frV>>rV[k];
	}
	fI.close();fQ.close();fV.close();frQ.close();frV.close();


    if (nth) {
      const string       dir_nth="/home/rgold/rt/newthermaltabs/";
      // const string       dir_nth="/home/rgold/codes/rt-git/newthermaltabs/";
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

    /************************* END OF LOOKUP TABLES *****************************/



	inited=true; //repeat once //RG:never repeated? if !(inited){}
    };
    
    if (false){ // CHECK lookup tables
    int inspect=0;
    int inspect_W=nWlen/2,inspect_T=0;
    int inspect_W_nth=nWlen_nth/2,inspect_T_nth=0;
    printf(YELLOW"[init.cpp]:"RESET" frequency=%E temperature=%E\n",12000.*pow(1.1, -nWlen/2+(inspect_W/2)),0.4*pow(1.05, inspect_T));
    printf(YELLOW"[init.cpp]:"RESET" frequency_nth=%E temperature_nth=%E\n",12000.*pow(1.25, -nWlen_nth/2+(inspect_W_nth/2)),0.4*pow(1.05, inspect_T_nth));
    printf(YELLOW"[init.cpp]:"RESET" exp(jI[%d][%d])=%E exp(jI_nth[%d][%d])=%E\n",inspect_T,inspect_W,exp(jI[inspect][inspect]),inspect_T_nth,inspect_W_nth,exp(jI_nth[inspect][inspect]));
    } // if () // CHECK lookup tables

}; // if(!inited){



/**********************************************
 * READ FLUID SIMULATION DUMP FILES FROM DISK *
 **********************************************/

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

 }; // for(i=-fdiff;i<=fdiff;i++) {




/**************************************************************
 * OVERWRITE DATA WITH AVERY'S FORCE-FREE TOY JET +RIAF MODEL *
 **************************************************************/
// int setup_averys_toyjet(doub rr, doub costh, float* uu[phlen][thlen][rlen][wdd], int isum);
 
if ( !strcmp(avery_toy_jet,"yes") ) { // consider use of PREPROCESSOR DIRECTIVE: #ifdef TOYJET ?

  // extern void setup_averys_toyjet(doub r, doub costh, doub rest[], doub a, 
  //                                 int isum, int r_idx, int th_idx, int phi_idx,
  //                                 doub gmunu[4][4], doub ginvmunu[4][4]);
  extern void setup_averys_toyjet(doub r, doub costh, doub rest[], doub a, 
                                  int isum, int r_idx, int th_idx, int phi_idx,
                                  doub gmunu[4][4], doub ginvmunu[4][4],
                                  doub FofPsi[][rlen], doub F[][rlen], doub &b_F_sq, doub &gamma);

  if (fdiff!=0) {
    fdiff=0;
    printf(YELLOW"[init.cpp]:"RESET"You called Avery's toyjet model with fdiff=%d. Will set fdiff=0 (model is stationary).\n",fdiff);
  }

  // #include "setup_averys_toyjet.cpp" //RG: can't define a fct inside another fct
  // prior to b_F_downr etc fix:
  // int rmin_debug=7; //0; 6: problem r~1.3 7: r~1.4 works
  int rmin_debug=0; //0; 25: problem r~1.3 26: r~2 works
  // int rmax_debug=rlen-1;//101; //rlen-1; //rlen; 262 works 263 problem: gamma nan
  int rmax_debug=ndd;//101; //rlen-1; //rlen; 262 works 263 problem: gamma nan
  int thmin_debug=0; //0;
  int thmax_debug=thlen-1; //thlen;
  int phmin_debug=0; //0;
  int phmax_debug=phlen; //phlen;
  printf(YELLOW"[init.cpp]: "RESET"smallest radius used in setup_averys_toyjet: %f\n",exp(coord[rmin_debug][0][0]));



  // ~> fct args
  doub r_fp=10;
  doub xi_toyjet=0.5;
  doub p_toyjet=2.-2.*xi_toyjet;
  printf(YELLOW"[init.cpp]: "RED"[HARDWIRE WARNING] r_fp=%f xi_toyjet=%f p_toyjet=%f defined in two places. Code assumes they are the same...\n"RESET,r_fp,xi_toyjet,p_toyjet);


  bool compute_F=false;
  if (compute_F) {
    
    //rmin_debug=50;
    // Build F(ψ) array first -> Needed for density in jet
    doub FofPsi[2][rlen];

    for(int x_index=rmin_debug;x_index<rmax_debug;x_index++) {
      doub x     = exp(coord[x_index][0][0]);
      doub z     = r_fp;
      doub costh = z/(x*x+z*z);
      doub cossq = costh*costh;
      doub sinsq = 1.-cossq;
      doub cos2th= cossq-sinsq;
      doub gtt_BL   = -(1.-2.*x/(x*x+asq*costh*costh)); // ~>BL[0][0]
      doub gtph_BL   = 2.*x/(x*x+asq*costh*costh) *a *sinsq; // ~>BL[0][3]
      doub gphph_BL   = (x*x +asq + 2.*x*asq/(x*x+asq*costh*costh)*sinsq)*sinsq;// ~>BL[0][0]
      // doub Sigma=x*x + asq*cossq; // Avoid confusion with Avery's sigma
      doub detg_BL=((x*x + asq*cossq)*(asq*(x - (x*x + asq*cossq)) + x*x*(2*x - (x*x + asq*cossq)) + asq*x*cos2th)*sinsq)/(asq - 2.*x + x*x);
      doub ginvtt_BL=-1./(x*x-2.*x+asq)*(x*x+asq+2.*x*asq*sinsq)/(x*x+asq*costh*costh);  // ~> iBL[0][0]
      doub ginvtph_BL=2.*x/(x*x+asq*costh*costh) *a * sinsq; // ~> iBL[0][3]
      doub sigma=1;
      doub OmegaOfPsi=pow(x,-3./2.);  OmegaOfPsi = max(OmegaOfPsi , pow(6.,-3/2.));   // floor@ISCO
      
      // BL metric: CHECK INVERSE: BL.iBL = 1 ?
      rhosq=x*x+asq*cossq;
      Delta=asq-2.*x+x*x;
      BL[0][0]=-(1. - 2.*x/rhosq);
      BL[0][1]=0.;
      BL[0][2]=0.;
      // BL[0][3]=Delta/rhosq * a *sinsq - sinsq/rhosq*(rsq+asq)*a; // does not look right...
      BL[0][3]=-2.*x/rhosq *a *sinsq;
      BL[1][0]=BL[0][1];
      BL[1][1]=rhosq/Delta;
      BL[1][2]=0.;
      BL[1][3]=0.;
      BL[2][0]=BL[0][2];
      BL[2][1]=BL[1][2];
      BL[2][2]=rhosq;
      BL[2][3]=0.;
      BL[3][0]=BL[0][3];
      BL[3][1]=BL[1][3];
      BL[3][2]=BL[2][3];
      // BL[3][3]=Delta/rhosq*a*sinsq + sinsq/rhosq*(rsq+asq); //does not look right...
      BL[3][3]=(x*x + asq + 2.*x*asq/rhosq*sinsq)*sinsq; // eq (3.11)
      
      // contravariant (upper indices) BL metric matrix
      // www.roma1.infn.it/teongrav/leonardo/bh/bhcap3.pdf
      iBL[0][0] = -(1./Delta) * ( x*x + asq + 2.*x*asq*sinsq/rhosq ); // eq (3.17)
      iBL[0][1] = 0.;
      iBL[0][2] = 0.;
      iBL[0][3] = -2.*x*a/Delta/rhosq; // eq (3.17)
      iBL[1][0] = iBL[0][1];
      iBL[1][1] = Delta/rhosq; // iBL[1][1] = 1 / BL[1][1]
      iBL[1][2] = 0.;
      iBL[1][3] = 0.;
      iBL[2][0] = iBL[0][2];
      iBL[2][1] = iBL[1][2];
      iBL[2][2] = 1./rhosq; // iBL[2][2] = 1/BL[2][2]
      iBL[2][3] = 0.;
      iBL[3][0] = iBL[0][3];
      iBL[3][1] = iBL[1][3];
      iBL[3][2] = iBL[2][3];
      iBL[3][3] = (Delta - asq*sinsq) / (rhosq*Delta*sinsq); // eq (3.17)

      

      /***********************
       * QUERY ANALYTIC MODEL see [setup_averys_toyjet.cpp] to compute F(Psi)
       ***********************/

      //printf(YELLOW"[init.cpp]: "RESET"Calling setup_averys_toyjet()() to get F. x=%f\n",x);

      // >~~~
      setup_averys_toyjet(x, costh, rest, a, x_index, x_index, 0, 0, BL, iBL, FofPsi, F, b_F_sq, gamma); // Only purpose here is to get F(Psi). 
      // <~~~

      //printf(YELLOW"[init.cpp]: "RESET"Exiting setup_averys_toyjet()(). x=%f\n",x);    
      
      //COPY ARRAY
      //doub F[2][rlen]; // Copy to be used for 2nd call to setup_averys_toyjet() ~> later FofPsi will be overwritten
      
      //FIXME: How to copy array?
      //std::copy ( FofPsi[2][rlen], FofPsi[2][rlen]+rlen, F[2][rlen].begin() );
      // simple cpp way:
      // std::vector<int> arr = {1, 2, 3, 4, 5};
      // std::vector<int> copy = arr;
      for (int i=0; i<rlen;i++) {
        F[0][i]=FofPsi[0][i];
        F[1][i]=FofPsi[1][i];
      }
      //COPY ARRAY

      //printf(YELLOW"[init.cpp]: "RESET"x= %f FofPsi[0][%d] = %f FofPsi[1][%d] = %f\n",x,x_index,FofPsi[0][x_index],x_index,FofPsi[1][x_index]);

    } //  for(int x_index=rmin_debug;x_index<rmax_debug;x_index++) {

  } // if (compute_F) {




  for(int theta_index=thmin_debug;theta_index<thmax_debug;theta_index++)
    for(int phi_index=phmin_debug;phi_index<phmax_debug;phi_index++)
      for(int r_index=rmin_debug;r_index<rmax_debug;r_index++){
  // for(int theta_index=0;theta_index<thlen-1;theta_index++)
  //   for(int phi_index=0;phi_index<phlen;phi_index++)
  //     for(int r_index=rmin_debug;r_index<rmax_debug;r_index++){

        int isum=r_index+theta_index+phi_index;

        //if (phi_index+theta_index==0) printf(YELLOW"[init.cpp]: "RESET"r_index=%d theta_index=%d phi_index=%d isum=%d\n",r_index,theta_index,phi_index,isum);

        doub r     = exp(coord[r_index][theta_index][0]);
        // doub r     = coord[r_index][theta_index][0]; 

        doub costh = coord[r_index][theta_index][1]; //RG:CHECK theta or cos(th)?
        doub z    = costh*r;
        doub Rc   = sqrt(r*r-z*z); //cylindrical radius
        doub Rb   = 20.;

        doub n0_th_RIAF  = 1.23e4; // 1.23x10^4 cm^-3 see main text after eq (4)
        doub n0_nth_RIAF = 3.8e2;  // 6.1 x10^2 cm^-3 see main text after eq (4)
        //doub n_th_RIAF  = n0_th_RIAF  * exp(-z*z/(2.*Rc*Rc)) * pow(r/Rb,-0.7);
        doub n_th_RIAF  = n0_th_RIAF * pow(r/Rb,-0.7);//  * exp(-z*z/(2.*rsq)) * pow(r/Rb,-0.7);
        doub n_nth_RIAF = n0_nth_RIAF * exp(-z*z/(2.*Rc*Rc)) * pow(r/Rb,-2.);

        if (isum==0) printf(YELLOW"[init.cpp]: "RED"n_th_RIAF=%e n0_th_RIAF=%e z=%e r=%e Rc=%e\n"RESET,n_th_RIAF,n0_th_RIAF,z,r,Rc);

        doub cossq = costh*costh;
        doub rsq = r*r;
        // a = 0.;
        doub asq = a*a;

        doub sinsq = 1.-cossq;
        doub sinth = sqrt(sinsq);

        // contravariant/inverse KS metric (see mathematica notebook)
        // CHECKED: THESE ARE INVERSE TO EACH OTHER TO MACHINE PRECISION
        KS[0][0] = -(r*(2.+r)+asq*cossq)/(rsq+asq*cossq);
        KS[0][1] = 2*r/(rsq+asq *cossq);
        KS[0][2] = 0; 
        KS[0][3] = 0;
        KS[1][0] = KS[0][1];
        KS[1][1] = 2.*(asq+(-2.+r)*r)/(asq+2.*rsq+asq*(cossq-sinsq)); // cos(2Theta)=cossq-sinsq // http://mathworld.wolfram.com/Double-AngleFormulas.html
        KS[1][2] = 0;
        KS[1][3] = a/(asq+rsq-asq*sinsq);
        KS[2][0] = KS[0][2];
        KS[2][1] = KS[1][2];
        KS[2][2] = 1/(rsq+asq*cossq);
        KS[2][3] = 0;
        KS[3][0] = KS[0][3];
        KS[3][1] = KS[1][3];
        KS[3][2] = KS[2][3];
        KS[3][3] = 1./sinsq/(asq+rsq-asq*sinsq);

        // covariant KS metric (see mathematica notebook)
        iKS[0][0] = -1. + 2.*r/(rsq+asq*cossq); 
        iKS[0][1] = 2.*r/(rsq+asq*cossq);
        iKS[0][2] = 0;
        iKS[0][3] = -2.*a*r*sinsq/(rsq+asq*cossq);
        iKS[1][0] = iKS[0][1];
        iKS[1][1] = 1. + 2.*r/(rsq+asq*cossq);
        iKS[1][2] = 0;
        iKS[1][3] = a*(-1. - 2.*r/(rsq+asq*cossq))*sinsq;
        iKS[2][0] = iKS[0][2];
        iKS[2][1] = iKS[1][2];
        iKS[2][2] = rsq+asq*cossq;
        iKS[2][3] = 0;
        iKS[3][0] = iKS[0][3];
        iKS[3][1] = iKS[1][3];
        iKS[3][2] = iKS[2][3];
        iKS[3][3] = sinsq*(asq+rsq+(2.*asq*r*sinsq)/(rsq+asq*cossq));


        //RG: USE rest array with setup_avery_toyjet
        //    Then set uu array inside init.cpp 
        //    (This avoids complications due to the weird uu-type)
        // printf(YELLOW"[init.cpp]: "RESET"Calling setup_averys_toyjet() again using F. r=%f costh=%f\n",r,costh);

        // >~~~
        //setup_averys_toyjet(r, costh, rest, a, isum, r_index, theta_index, phi_index, KS, iKS, FofPsi, F, b_F_sq, gamma); // query analytic model see [setup_averys_toyjet.cpp]
        //setup_averys_RIAF(r, costh, a);
        // <~~~



        // RIAF
        // (*uu[fdiff])[phi_index][theta_index][r_index][0] = n_th_RIAF*mp;
        // (*uu[fdiff])[phi_index][theta_index][r_index][1] = 8.1e9/mp;
        // (*uu[fdiff])[phi_index][theta_index][r_index][0] = n_th_RIAF;
        // (*uu[fdiff])[phi_index][theta_index][r_index][1] = 8.1e9*kb/mp;
        // (*uu[fdiff])[phi_index][theta_index][r_index][0] = n_th_RIAF;
        // (*uu[fdiff])[phi_index][theta_index][r_index][1] = 8.1e9*kb*n_th_RIAF*cc*cc;

        // P V = rho kT
        // P   = nkT
        // P   = (gamma-1) u
        // u   = P/(gamma-1) = nkT/(gamma-1)

        //printf(YELLOW"[init.cpp]: "RESET"n_th_RIAF=%e r=%f costh=%f\n",n_th_RIAF,r,costh);

        //        (*uu[fdiff])[phi_index][theta_index][r_index][0] = n_th_RIAF;
        // doub n_avery_RIAF=pow(r,-1.5);
        doub n_avery_RIAF=pow(0.5*r,-1.5);

        (*uu[fdiff])[phi_index][theta_index][r_index][0] = n_avery_RIAF;
        doub gamma=4./3.;
        (*uu[fdiff])[phi_index][theta_index][r_index][1] = n_avery_RIAF*8.1e9*kb/(gamma-1.); //RG:NOT NEEDED SEE rest[1] below

        (*uu[fdiff])[phi_index][theta_index][r_index][2] = 0.;
        (*uu[fdiff])[phi_index][theta_index][r_index][3] = 0.;
        // doub u_F_upt=;


        //define MKS to KS transformation for given coordinates
        for(i=0;i<4;i++)
          for(k=0;k<4;k++)
			MKStoKS[i][k]=dxdxp[r_ext_boundary_idx][nt][i][k];
        //RG: MKStoiKS

        //initialize with zeros metric tensor in MKS
        for(i=0;i<4;i++)
          for(k=0;k<4;k++)
			iMKS[i][k]=0.;
        //RG:2DO MKS[i][k]=0.;
        //compute metric tensor in MKS via a basic tensor transformation
        for(i=0;i<4;i++)
          for(k=0;k<4;k++)
			for(j=0;j<4;j++)
              for(m=0;m<4;m++)
                iMKS[i][k]+=MKStoKS[m][i]*iKS[m][j]*MKStoKS[j][k];
    //RG:2DO MKS[i][k]+=MKStoKS[m][i]*KS[m][j]*MKStoKS[j][k];



        // doub gtt=KS[0][0]; doub gtph=KS[1][3]; doub gphph=KS[3][3];
        doub gtt=iMKS[0][0]; doub gtph=iMKS[1][3]; doub gphph=iMKS[3][3];

        doub u_upphi = pow(Rc,-3./2.);   // Kepler's law (floor@ISCO)
        //u_upphi = max(u_upphi , pow(6.,-3/2.));   // floor@ISCO
        doub u_F_upt     = 1./sqrt( fabs(gtt)); //ZAMO observer

        (*uu[fdiff])[phi_index][theta_index][r_index][4] = u_F_upt; //WIP
        //(*uu[fdiff])[phi_index][theta_index][r_index][4] = 1; //WIP
        (*uu[fdiff])[phi_index][theta_index][r_index][5] = 0.;
        //(*uu[fdiff])[phi_index][theta_index][r_index][5] = -1e-2;
        (*uu[fdiff])[phi_index][theta_index][r_index][6] = 0.;
        (*uu[fdiff])[phi_index][theta_index][r_index][7] = 0.; // u_upphi;
        //(*uu[fdiff])[phi_index][theta_index][r_index][7] = 0.;
        (*uu[fdiff])[phi_index][theta_index][r_index][8] = 0.;
        (*uu[fdiff])[phi_index][theta_index][r_index][9] = 0.;
        doub plasma_beta=10.;
        doub B_sq       = 8.*PI/plasma_beta *n_th_RIAF *mp*cc*cc/6./r; // eq(3)
        (*uu[fdiff])[phi_index][theta_index][r_index][10]= 0.; //sqrt(B_sq);

        for(m=0;m<11;m++)
			rest[m]=(*uu[fdiff])[phi_index][theta_index][r_index][m];
			// rest[m]=(*uu[fdiff])[np][nt][r_ext_boundary_idx][m];

        if (isum==0) printf("rest[0]=%f rhonor=%f n_th_RIAF=%f\n",rest[0],rhonor,n_th_RIAF);


        rho=rest[0]*rhonor;                //physical density
        doub Te_RIAF=8.1e9;
        //doub Rb=20.;
        Te_RIAF *= pow(r/Rb,-0.84);
        rest[1] = Te_RIAF/mp/cc/cc*3*kb*rest[0];
        //Ttot=rest[1]*mp*cc*cc/3/kb/rest[0];//internal energy density
        u[0]=rest[4];                      //Lorentz factor
        u[1]=u[0]*rest[5];                 //u^r
        u[2]=u[0]*rest[6];                 //u^\theta
        u[3]=u[0]*rest[7];                 //u^\phi
        doub Bnor=sqrt(4.*PI*rhonor*mp)*cc;
        Bi[1]=Bnor*rest[8];                //magnetic field 3-vector
        Bi[2]=Bnor*rest[9];
        Bi[3]=Bnor*rest[10];

        // rho=rest[0]*rhonor;                //physical density
        // Ttot=rest[1]*mp*cc*cc/3/kb/rest[0];//internal energy density
        // u[0]=rest[4];                      //Lorentz factor
        // u[1]=u[0]*rest[5];                 //u^r
        // u[2]=u[0]*rest[6];                 //u^\theta
        // u[3]=u[0]*rest[7];                 //u^\phi
        // Bi[1]=Bnor*rest[8];                //magnetic field 3-vector
        // Bi[2]=Bnor*rest[9];
        // Bi[3]=Bnor*rest[10];


		//rho=rest[0]*rhonor;              //density
		//rest[1]*=mp*cc*cc/3/kb/(rest[0]);//temperature //RG:CHECK UNITS WHEN USING AVERY'S MODEL
		//rest[1]*=1./3./kb/(rest[0]);//temperature //RG:CHECK UNITS WHEN USING AVERY'S MODEL
		//rest[1]=8.1e9;//temperature //RG:CHECK UNITS WHEN USING AVERY'S MODEL
        // Ttot=8.1e9;

		// u[0]=rest[4];                    //4-velocity
		// u[1]=u[0]*rest[5];
		// u[2]=u[0]*rest[6];
		// u[3]=u[0]*rest[7];
        // //Bnor=sqrt(4.*PI*rhonor*mp)*cc;   //normalization of physical magnetic field, negative normalization is also allowed
		// Bi[1]=Bnor*rest[8];              //3-vector of magnetic field
		// Bi[2]=Bnor*rest[9];
		// Bi[3]=Bnor*rest[10];




        // Now that we know rest[] set/overwrite other variables used by astroray
		// for(int var=0;var<11;var++) {
        //   rest[var] = 0.; // OVERWRITE FLUID DATA WITH ZEROS
        //   //TEST
        //   if (isum==rmin_debug) printf(YELLOW"[init.cpp]: "RED"TEST\n"RESET);
        //   (*uu[fdiff])[phi_index][theta_index][r_index][var] = rest[var]; // OVERWRITE FLUID DATA (inverse of what's done elsewhere)
        //   if (rest[var]==0.) {
        //     //printf(YELLOW"[init.cpp]: "RED"rest[%d]=%f\n"RESET,var,rest[var]);
        //     //exit(1);
        //   }
        // }

		// rho=rest[0]*rhonor;              //density
        // // RG:FIXME? UNKNOWN TO THIS SCOPE...
		// // Ttot=rest[1]*mp*cc*cc/3/kb/rest[0];//temperature
		// u[0]=rest[4];                    //4-velocity
		// u[1]=u[0]*rest[5];
		// u[2]=u[0]*rest[6];
		// u[3]=u[0]*rest[7];

        // Bnor=sqrt(4.*PI*rhonor*mp)*cc; //normalization of physical magnetic field, negative normalization is also allowed

		// Bi[1]=Bnor*rest[8];              //3-vector of magnetic field
		// Bi[2]=Bnor*rest[9];
		// Bi[3]=Bnor*rest[10];

        // if (theta_index+phi_index==0 && r==50.) printf(YELLOW"[init.cpp]: "RESET"r=%e rho=%e u[0]=%e u[1]=%e Bi[1]=%e\n"RESET,r,rho,u[0],u[1],Bi[1]);

      } // r,theta,phi LOOP

  printf(YELLOW"[init.cpp]: "RESET"(*uu[0])[0][0][0][0]=%f\n",(*uu[fdiff])[0][0][0][0]);

  // for(int theta_index=thmin_debug;theta_index<thmax_debug;theta_index++) // only boundaries
  //   for(int phi_index=phmin_debug;phi_index<phmax_debug;phi_index++)
  //     for(int r_index=rmax_debug;r_index<rlen;r_index++){

  //   	//for(int var=0;var<11;var++)
  //   	for(int var=0;var<1;var++) //only rho
  //         (*uu[fdiff])[phi_index][theta_index][r_index][var] = 0.; // OVERWRITE FLUID DATA (inverse of what's done elsewhere)
  //       ;
  //     }

  if (count_nan_gamma!=0) {
    printf(YELLOW"[init.cpp]: "RED"THERE HAVE BEEN %d NANs in gamma\n"RESET,count_nan_gamma);
    printf(YELLOW"[init.cpp]: "RED"THE [SMALLEST,LARGEST] RADII WHERE NANs OCCURED ARE: [%f,%f]\n"RESET,smallest_radius_where_nan,largest_radius_where_nan);
  }
  if (count_cured_nan_gamma!=0) printf(YELLOW"[init.cpp]: "GREEN"%d NANs in gamma have been cured by different choice of gamma\n"RESET,count_cured_nan_gamma);

  // printf(YELLOW"[init.cpp]: "RED"EXITING\n"RESET); exit(1);

 } // AVERY's MODEL



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
printf(YELLOW"[init.cpp]: "RESET"Boundary between relaxed simulation data and radial extension zone: rcut=%f\n",rcut);
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
        // doub n_avery_RIAF=pow(0.5*r,-1.5);
        // (*uu[fdiff])[phi_index][theta_index][r_index][0] = n_avery_RIAF;
Bnor=sqrt(4.*PI*rhonor*mp)*cc; //normalization of physical magnetic field, negative normalization is also allowed

//extending radial temperature profile outwards as a powerlaw
for(k=ncut;k<ndd;k++){
	rtab[k]=rcut*pow((doub)rrmax/rcut,(k-ncut+(doub)1.)/(doub)(ndd-ncut));
	Tstab[k]=Ucon*pow((doub)rtab[k]/rcut,(doub)-Upo);
    //RG: ZERO OUT TEMPERATURES IN OUTER DOMAIN
    // if (k==ncut) printf(YELLOW"[init.cpp]: "RED"Zero out T in extended radial domain\n"RESET); 
    // Tstab[k]=0.;
};

int r_mdot_idx=60;                       //radial point where accretion rate is computed
doub rx=exp(coord[r_mdot_idx-1][0][0]),  //radius at r_mdot_idx
	 u0;                         //u^t
printf(YELLOW"[init.cpp]:"RED" HARDWIRED r[r_mdot_idx=%d]=%f (radius where mdot is evaluated)\n"RESET,r_mdot_idx,rx);
rate=0.;

//printf(YELLOW"[init.cpp]: "RESET"(*uu[0])[0][0][0][0]=%f (*uu[0])[0][0][10][0]=%f\n",(*uu[0])[0][0][0][0],(*uu[0])[0][0][10][0]);

//computing the accretion rate in code units
for(k=0;k<thlen-1;k++)
	for(i=0;i<phlen;i++)
      for(j=0;j<=2*fdiff;j++) { // RG: j -> t ?
			u0=(*uu[j])[i][k][r_mdot_idx-1][4];
			usp[j][i][k][0]=u0;
			usp[j][i][k][1]=u0*(*uu[j])[i][k][r_mdot_idx-1][5];//u^r
			usp[j][i][k][2]=u0*(*uu[j])[i][k][r_mdot_idx-1][6];//u^theta
			usp[j][i][k][3]=u0*(*uu[j])[i][k][r_mdot_idx-1][7];//u^phi
			uspKS[j][i][k]=0.;
			for(m=0;m<4;m++)
				uspKS[j][i][k]+=dxdxp[r_mdot_idx-1][k][1][m]*usp[j][i][k][m];//radial velocity in KS metric
			rate+=2./phlen*PI*(rx*rx+ asq*theta[r_mdot_idx-1][k]*theta[r_mdot_idx-1][k+1])*(*uu[j])[i][k][r_mdot_idx-1][0]*uspKS[j][i][k]*(theta[r_mdot_idx-1][k]-theta[r_mdot_idx-1][k+1]);

            if (isnan(rate)){ //RG: Check for NAN
              printf(YELLOW"[init.cpp]:"RESET"rate=%e,i=%d,j=%d,k=%d,rhonor=%e,rgrav=%e,fdiff=%d\n",rate,i,j,k,rhonor,rgrav,fdiff);
              printf(YELLOW"[init.cpp]:"RESET"rate=%e,i=%d,j=%d,k=%d,rx=%e,asq=%e,theta[r_mdot_idx-1][%d]=%e\n",rate,i,j,k,rx,asq,k+1,theta[r_mdot_idx-1][k+1]);
              printf(YELLOW"[init.cpp]:"RESET"rate=%e,uu[%d][%d][%d][r_mdot_idx-1][0]=%e,uspKS[j][i][k]=%e\n",rate,j,i,k,(*uu[j])[i][k][r_mdot_idx-1][0],uspKS[j][i][k]);
              printf(YELLOW"[init.cpp]:"RESET"rate=%e,dxdxp[r_mdot_idx-1][k][1][0]=%e,usp[j][i][k][0]=%e\n",rate,dxdxp[r_mdot_idx-1][k][1][0],usp[j][i][k][0]); 
              printf(YELLOW"[init.cpp]:"RESET"rate=%e,dxdxp[r_mdot_idx-1][k][1][1]=%e,usp[j][i][k][1]=%e\n",rate,dxdxp[r_mdot_idx-1][k][1][1],usp[j][i][k][1]);
              printf(YELLOW"[init.cpp]:"RESET"rate=%e,dxdxp[r_mdot_idx-1][k][1][2]=%e,usp[j][i][k][2]=%e\n",rate,dxdxp[r_mdot_idx-1][k][1][2],usp[j][i][k][2]); 
              printf(YELLOW"[init.cpp]:"RESET"rate=%e,dxdxp[r_mdot_idx-1][k][1][3]=%e,usp[j][i][k][3]=%e\n",rate,dxdxp[r_mdot_idx-1][k][1][3],usp[j][i][k][3]); 

              exit(1); /***************************************** EXIT *****************/
            }

      }; // for(k=0;k<thlen-1;k++) for(i=0;i<phlen;i++) for(j=0;j<=2*fdiff;j++) {


rate*=rhonor*rgrav*rgrav*cc*mp/(2*fdiff+1); // accretion rate in physical units

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

  // if (stNt==12107) printf(YELLOW"[init.cpp]: "RED"YO1: rz=%g rmin=%g stNt=%d\n"RESET,rz,rmin,stNt);

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

    //RG:FIXME DIFFERENT WHEN USING AVERY's MODEL... tp ->virial
    if ( !strcmp(avery_toy_jet,"yes") ) {
      // tp[stNt] = 8.1e9*pow(rr/20.,-0.84); 
      // te[stNt] = tp[stNt];
    }

    // printf(YELLOW"[init.cpp]: "RED"YO1: rz=%g rmin=%g status=%d stNt=%d\n"RESET,rz,rmin,status,stNt);

	if((rz<6) && (fl)) {    //catching electron temperature & Tp/Te ratio at 6M distance from the BH. Step size is so small that we don't iterate for 6M precisely
		fl=false;
		TpTe=tp[stNt]/te[stNt];
		Te6=te[stNt];
		printf(YELLOW"[init.cpp]:"RESET" fn=%d stN=%d r=%.2fM Tp/Te=%.2f Te=%.3e rate=%.3e he=%.3f rho=%.3e\n", fnum, stNt,rz, TpTe, Te6, rate*year/Msun,heat,rhonor);
	}

    // if (stNt==12107) printf(YELLOW"[init.cpp]: "RED"YO2: rz=%g rmin=%g status=%d stNt=%d\n"RESET,rz,rmin,status,stNt);


} // while (rz > 1.001*rmin) {

maxT=ts[stNt];minT=ts[0];

if(stNt>maxst){            //check if solver exceeded the number of steps
	printf(YELLOW"[init.cpp]:"RED"Temperature solver error. Exceeded maximum steps allowed! Increase maxst in [global_variables.cpp]\nEXITING\n"RESET);
	exit(-1);
};
gsl_odeiv_evolve_free (ez);//free memory
gsl_odeiv_control_free (cz);
gsl_odeiv_step_free (sz);



// RG: CODE EXITS WHEN UNCOMMENTING BELOW CODE BLOCK
// if (turn_off_radial_extension=true) {
//   t_init += (clock() - t_b4_init) / (doub)CLOCKS_PER_SEC;
  
//   printf("\n"YELLOW"[init.cpp]:"RESET" Exiting [init.cpp]\n");

//   return 0; // init()
//  }



/*************************
 * OUTER RADIAL BOUNDARY *
 *************************/

printf(YELLOW"[init.cpp]: "RED"SETTING OUTER RADIAL EXTENSION...\n"RESET);

//compute quantities in the co-moving frame at the convergence boundary (given by rr=rcut=const)
r_ext_boundary_idx=ncut-1; //point just inside convergence boundary
rr=rcut;    

/****************************************************************************/
for(nt=0;nt<thlen;nt++){ //for all POLAR angles at the convergence boundary *
	costh=theta[r_ext_boundary_idx][nt];
	cossq=costh*costh;
	sinsq=1-cossq;
    sinth=sqrt(sinsq);
	rsq=rr*rr;
	rhosq=rsq+asq*cossq;
	Delta=rsq-2.*rr+asq; // Lambda_squared in http://grwiki.physics.ncsu.edu/wiki/Kerr_Black_Hole
	doub temp=2.*rr/rhosq;

	//-,+,+,+ signature covariant KS metric matrix
    //compare to mathematica notebook... why iKS? This is the covariant metric
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
			MKStoKS[i][k]=dxdxp[r_ext_boundary_idx][nt][i][k];
    //RG: MKStoiKS

	//initialize with zeros metric tensor in MKS
	for(i=0;i<4;i++)
		for(k=0;k<4;k++)
			iMKS[i][k]=0.;
    //RG:2DO MKS[i][k]=0.;
	//compute metric tensor in MKS via a basic tensor transformation
	for(i=0;i<4;i++)
		for(k=0;k<4;k++)
			for(j=0;j<4;j++)
				for(m=0;m<4;m++)
					iMKS[i][k]+=MKStoKS[m][i]*iKS[m][j]*MKStoKS[j][k];
    //RG:2DO MKS[i][k]+=MKStoKS[m][i]*KS[m][j]*MKStoKS[j][k];
	
    /****************************************************/
	for(np=0;np<phlen;np++){ //for all AZIMUTHAL angles *

		//read fluid quantities at a given point
        // DON'T DO THIS WHEN USING ANALYTIC (e.g. AVERY's) MODEL...
        // if ( strcmp(avery_toy_jet,"yes") ) { // Consider use of PREPROCESSOR DIRECTIVE: #ifdef TOYJET ?
        for(m=0;m<11;m++)
			rest[m]=(*uu[fdiff])[np][nt][r_ext_boundary_idx][m];

        //RG:FIXME
        // printf(YELLOW"[init.cpp]: "RED"ZERO OUT RHO IN OUTER BOUNDARY\n"RESET);
        // rest[0]=0.; // TURN OFF (FOR AVERY'S MODEL)
        
        // } // if ( strcmp(avery_toy_jet,"yes") ) { 

		rho=rest[0]*rhonor;              //density
		rest[1]*=mp*cc*cc/3/kb/(rest[0]);//temperature //RG:CHECK UNITS WHEN USING AVERY'S MODEL
		u[0]=rest[4];                    //4-velocity
		u[1]=u[0]*rest[5];
		u[2]=u[0]*rest[6];
		u[3]=u[0]*rest[7];
		Bi[1]=Bnor*rest[8];              //3-vector of magnetic field
		Bi[2]=Bnor*rest[9];
		Bi[3]=Bnor*rest[10];
        
        // AVERY's model
        if (np+nt==0) printf(YELLOW"[init.cpp]: "RESET"Setting AVERY's model in radial extension zone...\n");

        // doub n_avery_RIAF=pow(0.5*exp(rcut),-1.5);
        doub n_avery_RIAF=pow(0.5*rcut,-1.5);
        // doub n_avery_RIAF=pow(0.5*coord[ncut][0][0],-1.5);
        rho=n_avery_RIAF;

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



// for(int theta_index=0;theta_index<thlen;theta_index++) // only boundaries
//   for(int phi_index=0;phi_index<phlen;phi_index++)
//     for(int r_index=141;r_index<rlen;r_index++){
      
//       //for(int var=0;var<11;var++)
//       for(int var=0;var<1;var++) //only rho
//         (*uu[fdiff])[phi_index][theta_index][r_index][var] = 0.; // OVERWRITE FLUID DATA (inverse of what's done elsewhere)
//       ;
//     }


t_init += (clock() - t_b4_init) / (doub)CLOCKS_PER_SEC;

printf("\n"YELLOW"[init.cpp]:"RESET" Exiting [init.cpp]\n");

return 0; // init()
}


