// void setup_averys_toyjet(doub r, doub costh, doub rest[], doub a, 
//                          int isum, int r_idx, int th_idx, int phi_idx,
//                          doub gmunu[4][4], doub ginvmunu[4][4]) {
void setup_averys_toyjet(doub r, doub costh, doub rest[], doub a, 
                         int isum, int r_idx, int th_idx, int phi_idx,
                         doub gmunu[4][4], doub ginvmunu[4][4], doub FofPsi[2][rlen], doub F[2][rlen], doub &b_F_sq, doub &gamma) {

  /************************************************************************************
   * Set input data for ASTRORAY (overwrite simulation data) 
   * to toy-jet+RIAF model by Broderick & Loeb 2009 
   * The Astrophysical Journal, 697:1164–1179, 2009 June 1
   * "IMAGING THE BLACK HOLE SILHOUETTE OF M87: 
   *  IMPLICATIONS FOR JET FORMATION AND BLACK HOLE SPIN"
   * http://dx.doi.org/10.1088/0004-637X/697/2/1164
   * http://adsabs.harvard.edu/abs/2009ApJ...697.1164B 
   * http://iopscience.iop.org/0004-637X/697/2/1164/ <- Reference for equation nrs
   * http://arxiv.org/abs/0812.0366
   * http://inspirehep.net/record/804072
   ************************************************************************************/

  bool JET  = false; // INCLUDE JET?
  bool RIAF = false; // INCLUDE RIAF? //WHATS WRONG WITH RIAF GET NANs in I -> rest[]=0


// single fluid simulation dump files contains: 
// uu[?]:     [0    1   2     3                4    5    6        7      8    9        10   ]
//             rho, u, -u^t, -T^r_t/(rho u^r), u^t, v^r, v^theta, v^phi, B^r, B^theta, B^phi

  // for(int theta_index=0;theta_index<thlen-1;theta_index++) // only boundaries
  //   for(int phi_index=0;phi_index<phlen;phi_index++)
  //     for(int r_index=0;r_index<rlen;r_index++){

  //       int isum=r_index+theta_index+phi_index;

  // if (isum==0) printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET" Get coordinates in here...\n");
  // doub r     = exp(coord[r_index][theta_index][0]);
  // doub r     = coord[r_index][theta_index][0]; 
  // doub costh = coord[r_index][theta_index][1]; //RG:CHECK theta or cos(th)?
  // rr=r;

  if (isum==0) cout <<YELLOW"[setup_averys_toyjet.cpp]:"RESET" WIP What about dPhi here...?!"<<endl;
  // doub phi   = ((doub) phi_index)/phlen / 2./PI; //RG: dPhi OFFSET?

  for (int mu=0; mu<=3;mu++) for (int nu=0; nu<=3;nu++) {
      //printf(YELLOW"[setup_averys_toyjet.cpp]: "RESET"[%d][%d]: gmunu=%f ginvmunu=%f\n",mu,nu,gmunu[mu][nu],ginvmunu[mu][nu]);
    }

  // costh=theta[r_index][theta_index];
  doub cossq=costh*costh;
  doub sinsq=1-cossq;
  doub sinth=sqrt(sinsq); //how accurate is that?
  doub rsq=r*r;
  doub rhosq=rsq+asq*cossq; // akas Sigma
  doub Delta=rsq-2.*r+asq; // Lambda_squared in http://grwiki.physics.ncsu.edu/wiki/Kerr_Black_Hole
  doub temp=2.*r/rhosq;
  doub z=r*costh;
  doub Rc=sqrt(rsq-z*z); //cylindrical radius

  // METRIC

  // // Boyer-Lindquist metric
  doub BL[4][4]; doub iBL[4][4];
  BL_metric(BL, iBL, r, costh); // SET METRIC AND ITS INVERSE TO BOYER-LINDQUIST METRIC

  // SET THE METRIC AND ITS INVERSE
  // KS_metric(gmunu, ginvmunu, r, costh); // SET METRIC AND ITS INVERSE TO KERR-SCHILD METRIC


  // if (isum==50) printf(YELLOW"[setup_averys_toyjet.cpp]: KS[0][0]=%e iKS[0][0]=%e gmunu[0][0]=%e ginvmunu[0][0]=%e\n"RESET,KS[0][0],iKS[0][0],gmunu[0][0],ginvmunu[0][0]);

  double kronecker[4][4];
  // int kronecker[4][4];
  for (int i=0; i<=3; i++) for (int j=0; j<=3; j++) { if (i==j) {kronecker[i][j]=1;} else kronecker[i][j]=0;}

  // CHECK METRIC (INVERSE)
  if (isum==100) {
    doub delta[4][4]; // if everything is right ~> identity matrix
    //int delta[4][4]; // if everything is right ~> identity matrix

    for (int mu=0; mu<=3; mu++) for (int nu=0; nu<=3; nu++) {
        delta[mu][nu]=0.;
      }
    for (int mu=0; mu<=3; mu++) for (int nu=0; nu<=3; nu++) for (int kappa=0; kappa<=3; kappa++) {
          //            for (int nu=0; mu<=3; mu++ and nu++) {
          delta[mu][nu] += BL[mu][nu]*iBL[mu][nu]; // delta *IS* diagonal but closer to 4*kronecker :-S
          // delta[mu][nu] += KS[mu][kappa]*iKS[kappa][nu]; 
          //delta[mu][nu] += gmunu[mu][kappa]*ginvmunu[kappa][nu];
          // delta[mu][nu] += kronecker[mu][kappa]*kronecker[kappa][nu];
        }
    
    for (int mu=0; mu<=3; mu++) for (int nu=0; nu<=3; nu++) { 
        //printf("delta[%d][%d]=%e\n",mu,nu,delta[mu][nu]);
        if (fabs(delta[mu][nu]-kronecker[mu][nu])>1e-10) printf("delta[%d][%d]=%e\n",mu,nu,delta[mu][nu]);
      }

  } // if (isum==...) {
  
  
  /********************************************************************************************/
  // float *rho    = (uu[0])[phi_index][theta_index][r_index][0]; // rest-mass density
  // float *u      = (uu[0])[phi_index][theta_index][r_index][1]; // internal energy
  // // float *rho = (uu[0])[phi_index][theta_index][r_index][2];
  // // float *rho = (uu[0])[phi_index][theta_index][r_index][3]; // RG: not required by ASTRORAY I believe...
  // float *u0     = (uu[0])[phi_index][theta_index][r_index][4]; // u^t? // u0 -> u_0 -> W
  // float *vr     = (uu[0])[phi_index][theta_index][r_index][5]; // Avery gives ur=u0*vr
  // float *vtheta = (uu[0])[phi_index][theta_index][r_index][6]; // Avery gives utheta=u0*vtheta
  // float *vphi   = (uu[0])[phi_index][theta_index][r_index][7]; // Avery gives uphi=u0*vphi
  // float *Br     = (uu[0])[phi_index][theta_index][r_index][8];
  // float *Btheta = (uu[0])[phi_index][theta_index][r_index][9];
  // float *Bphi   = (uu[0])[phi_index][theta_index][r_index][10];
  /********************************************************************************************/
  
  
  /*********************************
   * FREE PARAMETERS / USER CHOICE *
  /*********************************/
  doub n0_th_RIAF  = 1.23e4; // 1.23x10^4 cm^-3 see main text after eq (4)
  doub n0_nth_RIAF = 3.8e2;  // 6.1 x10^2 cm^-3 see main text after eq (4)
  doub Te_RIAF=8.1e9;        // Te_0 = 8.1x10^9K in eq (2) see main text after eq (4) actually not used here but in [evalpointzero.cpp]
  doub Rb=20; // Needed to reproduce for M87 (disk subdominant@7mm), see main text after eq (4)
  doub plasma_beta = 10;     // beta=10 see main text between eq (3) and (4)
  doub xi_toyjet = 0.5; // eq.(5) xi=0: cylindrical xi=1: conical determines collimation rate
  doub r_fp = 10.; // foot point of the jet where particles are injected
  doub p_toyjet = 2.-2.*xi_toyjet; // powerlaw-index in stream fct eq (A12) determines collimation rate 

  doub psi_toyjet = pow(r,p_toyjet) * (1. - costh); // stream fct eq.(5),(A42)
  // Omega(Psi): See main text after eq (5):
  // "The structure of the magnetic field and the velocity
  // of the outflow depend upon the angular velocity at the
  // field-footprint, Ω(ψ), which is a function solely of ψ. "

  //RG:HOWEVER BETWEEN eq(6) and eq(7):
  // "We define Ω(ψ) in the equatorial plane, choose it to be the Keplerian angular velocity out-
  // side of the ISCO, and fix it to the ISCO value inside. This
  // boundary is usually what demarcates the “jet”, separat-
  // ing it from the surrounding disk wind."

  doub OmegaOfPsi = pow(r,-3./2.);   // Kepler's law (floor@ISCO)
  OmegaOfPsi = max(OmegaOfPsi , pow(6.,-3/2.));   // floor@ISCO

  doub sigma=-1.; if (isum==0)printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET"[CHECK:]sigma=%f\n",sigma); // sigma=u_F.u_F=+/-1 see between eq.(6),(7)


  /**********
   * METRIC *
   **********/

  // covariant metric components
  //RG: USE METRIC & ITS INVERSE AS PASSED TO FCT
  doub gtt=gmunu[0][0];doub gtph=gmunu[0][3];doub gphph=gmunu[3][3];
  doub gtr=gmunu[0][1];doub grth=gmunu[1][2];doub grph=gmunu[1][3];
  doub gtth=gmunu[0][2]; doub gthph=gmunu[2][3];
  doub grr=gmunu[1][1]; doub gthth=gmunu[2][2];

  // doub detg=-sinsq*pow(rsq+asq*cossq,2); // See e.g. arXiv:0706.0622v3 eq(70)
  doub detg_BL = BL[0][0] * (BL[1][1]*BL[2][2]*BL[3][3]) + BL[0][3] * (BL[1][1]*BL[2][2]*BL[0][3]);
  //RG:inverse
  //doub detginv_KS =  -4./sinsq/pow(asq+2.*rsq+asq*(cossq-sinsq),2); //det(ginv) // cos(2Theta)=cossq-sinsq 
  doub detg_KS = (rsq + asq*cossq)*(-asq*sinsq - rsq*sinsq + asq*sinsq*sinsq); // DET OF COVARIANT KS METRIC -> MATHEMATICA

  doub detg=detg_KS;

  // FLAT METRIC
  // doub gtt=-1.;doub gtph=0.;doub gphph=rsq*sinsq;
  // doub grr=1.;doub gthth=rsq;
  // doub detg=1.;

  // contravariant metric components
  doub ginvtt   = ginvmunu[0][0];
  doub ginvtr   = ginvmunu[0][1];
  doub ginvrr   = ginvmunu[1][1]; // g^rr   = 1/g_rr (in BL)
  doub ginvthth = ginvmunu[2][2]; // g^thth = 1/g_thth (in BL)
  doub ginvphph = ginvmunu[3][3];
  doub ginvtph  = ginvmunu[0][3];
  doub ginvrph  = ginvmunu[1][3];

  // JET solution 
  // $ u^t_F $ eq.(A43) Does not assume BL. Generally covariant.
  doub u_F_upt     = 1./sqrt( fabs(gtt + 2.*gtph*OmegaOfPsi + gphph*pow(OmegaOfPsi,2)) ); // eq.(A43)

  // BUT u_F.u_F=sigma -> u_F^t = Omega gtph/gtt +/- sqrt( (Omega gtph/gtt)^2 -Omega^2 gphph/gtt )
  doub u_F_upr     = 0.;                 //eq (6)
  doub u_F_uptheta = 0.;                 //eq (6)
  doub u_F_upphi   = u_F_upt*OmegaOfPsi; // follows from Omega definition just before eq.(A30)

  doub u_F_squared=sigma; //after eq.(A7)

  // $ b^\mu_F $ eq.(8,A44)
  doub b_F_downphi = -2. * OmegaOfPsi * psi_toyjet * u_F_upt; // covariant
  doub b_F_downt   = -b_F_downphi * OmegaOfPsi; // covariant

  // $ b^\mu_F $ eq.(7,A44) and some simple algebra
  doub b_F_upr     = -sigma * pow(r,p_toyjet)*sinth / (u_F_upt * sqrt(-detg));
  doub b_F_uptheta = -sigma*p_toyjet * pow(r,p_toyjet-1)*(1.-costh) / (u_F_upt * sqrt(-detg));

  // doub b_F_downr     = grr*b_F_upr; //RG:FIXME: ASSUMES BL! GENERALIZE!
  // doub b_F_downr     = gtr*b_F_upt + grr*b_F_upr + grth*b_F_uptheta + grph*b_F_upphi;
  //RG:CHECK covariant?
  doub b_F_downr     = gtr*( ginvtt*b_F_downt + ginvtph*b_F_downphi + ginvtr*b_F_downr ) + grr*b_F_upr + grph*( ginvtph*b_F_downt + ginvrph*b_F_downr + ginvphph*b_F_downphi);

  //RG: Not covariant
  doub b_F_upphi   = ginvtph * b_F_downt + ginvphph * b_F_downphi + b_F_downr*ginvrph; // last term vanishes for BL
  doub b_F_upt     = ginvtt * b_F_downt + ginvtph * b_F_downphi + ginvtr*b_F_downr; // last term vanishes for BL

  doub b_F_downtheta = gtth*b_F_upt + grth*b_F_upr + gthth*b_F_uptheta + gthph*b_F_upphi; // All but 3rd term vanish for BL&KS

  // b_F^2 mostly cancels with beta^2 except for density... //RG:BL->KS ok
  //doub b_F_sq = b_F_downr*b_F_upr + b_F_downtheta*b_F_uptheta + b_F_downphi/u_F_upt*(u_F_upt*b_F_upphi - u_F_upphi*b_F_upt); // eq(A40) ~~~> [double-checked]
  // Now passed by reference via fct arg
  b_F_sq = b_F_downr*b_F_upr + b_F_downtheta*b_F_uptheta + b_F_downphi/u_F_upt*(u_F_upt*b_F_upphi - u_F_upphi*b_F_upt); // eq(A40) ~~~> [double-checked]

  doub b_F_sq_v2 = b_F_downt*b_F_upt + b_F_downr*b_F_upr + b_F_downtheta*b_F_uptheta + b_F_downphi*b_F_upphi; // eq(A40) ~~~> [double-checked] agrees with b_F_sq to within truncation error
 
  if ( fabs(1.-b_F_sq/b_F_sq_v2) > 1e-10 ) printf(YELLOW"[setup_averys_toyjet.cpp]: "RED"Something looks inconsistent with b_F_sq ...\n"RESET);

  doub beta = sigma * b_F_upt / (b_F_sq * u_F_upt); // eq.(10,A46) ~~~> [double-checked]


  // sometimes: beta*beta*b_Fsq >> sigma  => gamma=nan !
  /************************* GAMMA (NAN WATCH) **********************/
  //doub gamma = -sigma / sqrt( -(sigma + beta*beta*b_F_sq ) ); // eq.(A24) ~~~> [double-checked]
  // Now passed by reference via fct arg
  gamma = -sigma / sqrt( -(sigma + beta*beta*b_F_sq ) ); // eq.(A24) ~~~> [double-checked]
  /******************************************************************/


  if isnan(gamma) { //FIXME: need to recompute other things too
      // printf(YELLOW"[setup_averys_toyjet.cpp]: "RED"gamma=nan ! -> change sigma ...WIP...\n"RESET);
      sigma         *= -1;
      beta          *= -1;
      b_F_upr       *= -1;
      b_F_uptheta   *= -1;
      b_F_downr     *= -1;
      b_F_downtheta *= -1;
      // b_F_sq invariant under sigma*=-1 [CHECK]
      b_F_sq = b_F_downr*b_F_upr + b_F_downtheta*b_F_uptheta + b_F_downphi/u_F_upt*(u_F_upt*b_F_upphi - u_F_upphi*b_F_upt); // eq(A40) ~~~> [double-checked]
      
      /************************* GAMMA (nan watch) **********************/
      gamma = -sigma / sqrt( -(sigma + beta*beta*b_F_sq ) ); // eq.(A24) ~~~> [double-checked]
      /******************************************************************/

      if isnan(gamma) {
          printf(YELLOW"[setup_averys_toyjet.cpp]: "RED"gamma=nan for both choices of sigma r=%f cos(theta)=%f\n...BETTER DIE...\n"RESET,r,costh);
          printf(YELLOW"[setup_averys_toyjet.cpp]: "RED"beta=%f b_F_sq=%f beta^2*b_F^2=%f!\n"RESET,beta,b_F_sq,beta*beta*b_F_sq);
          printf(YELLOW"[setup_averys_toyjet.cpp]: "RED"gtt=%f gtph=%f gphph=%f OmegaOfPsi=%f!\n"RESET,gtt,gtph,gphph,OmegaOfPsi);
          //doub u_F_upt     = 1./sqrt( fabs(gtt + 2.*gtph*OmegaOfPsi + gphph*pow(OmegaOfPsi,2)) ); // eq.(A43)

          // printf("b_F_sq=%f b_F_sq_v2=%f\n",b_F_sq,b_F_sq_v2); //LOOKS OK

          //exit(1);
          count_nan_gamma++;
          if (r < smallest_radius_where_nan) smallest_radius_where_nan = r;
          if (r > largest_radius_where_nan)  largest_radius_where_nan  = r;
          // printf(YELLOW"[setup_averys_toyjet.cpp]: "RED"gamma=nan for both choices of sigma r=%f cos(theta)=%f!\n"RESET,r,costh);
          // printf(YELLOW"[setup_averys_toyjet.cpp]: "RED"beta=%f b_F_sq=%f beta^2*b_F^2=%f!\n"RESET,beta,b_F_sq,beta*beta*b_F_sq);

          // ZERO OUT ALL VARIABLES
          // for (int var=0; var<=10; var++) rest[var] = 0.;
        }
      else {
        count_cured_nan_gamma++;
        // printf(YELLOW"[setup_averys_toyjet.cpp]: "GREEN"gamma cured by choice of sigma=%f at r=%f costh=%f!\n"RESET,sigma,r,costh);
        // printf(YELLOW"[setup_averys_toyjet.cpp]: "RESET"beta=%f b_F_sq=%f beta^2*b_F^2=%f!\n",beta,b_F_sq,beta*beta*b_F_sq);
      }

    } // if isnan(gamma) {

  // JET: u_F,b_F -> u,b : eqs (9,11,A46)
  doub u_RIAFtoyjet[4];
  u_RIAFtoyjet[0] = gamma * (u_F_upt     + beta*b_F_upt);
  u_RIAFtoyjet[1] = gamma * (u_F_upr     + beta*b_F_upr);
  u_RIAFtoyjet[2] = gamma * (u_F_uptheta + beta*b_F_uptheta);
  u_RIAFtoyjet[3] = gamma * (u_F_upphi   + beta*b_F_upphi);

  // DISK
  // doub z=r*costh; //declared earlier

  // eqs.(1,4)
  doub n_th_RIAF  = n0_th_RIAF  * exp(-z*z/(2.*Rc*Rc));
  doub n_nth_RIAF = n0_nth_RIAF * exp(-z*z/(2.*Rc*Rc));

  if ( r >= Rb ) {
    n_th_RIAF  *= pow(r/Rb,-0.7);
    n_nth_RIAF *= pow(r/Rb,-2);  
  }

  //doub B_sq       = 8.*PI/plasma_beta *n_th_RIAF *   cc*cc/6./r; // eq(3) 
  doub B_sq       = 8.*PI/plasma_beta *n_th_RIAF *mp*cc*cc/6./r; // eq(3) 

  doub B_RIAF_phi = sqrt(B_sq); // purely toroidal

  // number density RIAF (DISK)
  //doub rho_RIAF = (n_th_RIAF + n_nth_RIAF)*mp;
  doub rho_RIAF = (n_th_RIAF + n_nth_RIAF);
  //doub rho_RIAF = (n_th_RIAF + n_nth_RIAF);//*mp;

  // JET (Sec. ?) & Appendix A <-> free fct F  eq.(13)
  //RG: DO THIS HERE NOT OUTSIDE setup_averys_toyjet()



  if (true) { //RG:FIXME NEED TO EVALUATE QUANTITIES IN DIFFERENT LOCATIONS ~> Can't do it here easily... ALTERNATIVE: RETURN b_F_sq,gamma

    // { //COPY ARRAY
    // //doub F[2][rlen]; // Copy to be used for 2nd call to setup_averys_toyjet() ~> later FofPsi will be overwritten

    // //FIXME: How to copy array?
    // //std::copy ( FofPsi[2][rlen], FofPsi[2][rlen]+rlen, F[2][rlen].begin() );
    // // simple cpp way:
    // // std::vector<int> arr = {1, 2, 3, 4, 5};
    // // std::vector<int> copy = arr;
    // for (int i=0; i<rlen;i++) {
    //   F[0][i]=FofPsi[0][i];
    //   F[1][i]=FofPsi[1][i];
    // }
    // } //COPY ARRAY

  // Build F(ψ) array first -> Needed for density in jet
  //doub FofPsi[2][rlen];

  //for(int x_index=0;x_index<rlen;x_index++) {
    
    int x_index=r_idx;
    doub xi_toyjet=0.5;
    doub p_toyjet=2.-2.*xi_toyjet;
    doub x_fp     = exp(coord[x_index][0][0]);
    doub z_fp     = r_fp; // spherical <-> cylindrical
    doub costh_fp = z_fp/sqrt(x_fp*x_fp + z_fp*z_fp);

    doub ntilde = (-(x_fp*x_fp+z_fp*z_fp)*sinsq/2/r_fp/r_fp); // between eqs (12,13) Broderick & Loeb 2009 

    FofPsi[0][x_index]  = pow(x_fp,p_toyjet) * (1. - costh_fp); // Set ψ
    FofPsi[1][x_index]  = gamma*ntilde/u_F_upt/b_F_sq*(iBL[0][0]+BL[0][3]/OmegaOfPsi); // Set F(ψ) according to eq(13) Broderick & Loeb 2009

    // if (x_index==100) {
    //   printf(YELLOW"[setup_averys_toyjet.cpp]: "RESET"FofPsi[0][%d] = %f FofPsi[1][%d] = %f\n",x_index,FofPsi[0][x_index],x_index,FofPsi[1][x_index]);
    // printf(YELLOW"[setup_averys_toyjet.cpp]: "RESET"gamma=%f ntilde=%f u_F_upt=%f b_F_sq=%f iBL[0][0]=%f BL[0][3]=%f OmegaOfPsi\n=%f\n",gamma,ntilde,u_F_upt,b_F_sq,iBL[0][0],BL[0][3],OmegaOfPsi);
    // // ~> [init.cpp]: FofPsi[0][0]=-8.832118 FofPsi[1][0]=nan
    // }
    
    //}  // for(int x_index=0;x_index<rlen-1;x_index++) { 

  }  // if (false/true) {


  doub F_intp = F_intp_lin(FofPsi, psi_toyjet); // INTERPOLATE F(Psi)
  //doub F=1.;

  // WITHOUT THIS LINE WE GET SEG-FAULT AFTER EXITING [init.cpp]
  for (int var=0; var<=10; var++) rest[var] = 0.; //AVOID UNINITIALIZED VALUES?!

  if (RIAF) {
    rest[0] += rho_RIAF; // careful with rest[0]=0 -> rest[1]
    //if (rest[0]==0.) printf(YELLOW"[init.cpp]: "RESET"RIAF model selected rest[0]=%f r=%f costh=%f\n",rest[0],r,costh);
  }
  if (JET) {
    if (isum==0) printf(YELLOW"[init.cpp]"RESET"JET model selected\n");
    rest[0] += F_intp * fabs( u_F_upt * b_F_sq / gamma / (ginvtt + ginvtph/OmegaOfPsi) * (1- exp(-pow(r/r_fp,2)/2.)) ); //eq.(13) F=1?
  }

  //rest[0]  = rho_RIAF; // debug: RIAF only

  // float *u0     = (uu[0])[phi_index][theta_index][r_index][4]; // u^t? // u0 -> u_0 -> W
  // float *vr     = (uu[0])[phi_index][theta_index][r_index][5]; // Avery gives ur=u0*vr
  // float *vtheta = (uu[0])[phi_index][theta_index][r_index][6]; // Avery gives utheta=u0*vtheta
  // float *vphi   = (uu[0])[phi_index][theta_index][r_index][7]; // Avery gives uphi=u0*vphi

  //RG:FIXME
  //for(int u_idx=4;u_idx<=7;u_idx++) (*uu[0])[phi_index][theta_index][r_index][4] = u_RIAFtoyjet[u_idx];
  
  rest[1] = Te_RIAF/(mp*cc*cc/3/kb/rest[0]); //Ttot=rest[1]*mp*cc*cc/3/kb/rest[0];//internal energy density

  if (RIAF) {
    rest[4] += 1.; // WIP
    rest[5] += 0.;
    rest[6] += 0.;
    rest[7] += pow(max(r,6.),-0.5); // KEPLERIAN
    // B
    rest[8] = 0.;
    rest[9] = 0.;
    rest[10]= B_RIAF_phi;
    // TEST
    // if (isum==0) cout<<"TESTING!"<<endl;
    // rest[8]  = 0.;
    // rest[9]  = 0.;
    // rest[10] = 0.;
    // rest[8]  = gamma * (b_F_upr - sigma*b_F_sq*beta*u_F_upr);
    // rest[9]  = gamma * (b_F_uptheta - sigma*b_F_sq*beta*u_F_uptheta);
    // rest[10] = gamma * (b_F_upphi - sigma*b_F_sq*beta*u_F_upphi);
  }

  if (JET) {
    rest[4] += gamma*(u_F_upt + beta*b_F_upt);
    rest[5] += gamma*(u_F_upr + beta*b_F_upr);
    rest[6] += gamma*(u_F_uptheta + beta*b_F_uptheta);
    rest[7] += gamma*(u_F_upphi + beta*b_F_upphi);
    // Magnetic field eq (11,A47)
    rest[8]  = gamma * (b_F_upr - sigma*b_F_sq*beta*u_F_upr);
    rest[9]  = gamma * (b_F_uptheta - sigma*b_F_sq*beta*u_F_uptheta);
    rest[10] = gamma * (b_F_upphi - sigma*b_F_sq*beta*u_F_upphi);
  }

  // (*uu[0])[phi_index][theta_index][r_index][8]  = gamma * (b_F_upr - sigma*b_F_sq*beta*u_F_upr);
  // (*uu[0])[phi_index][theta_index][r_index][9]  = gamma * (b_F_uptheta - sigma*b_F_sq*beta*u_F_uptheta);
  // (*uu[0])[phi_index][theta_index][r_index][10] = gamma * (b_F_upphi - sigma*b_F_sq*beta*u_F_upphi) + B_RIAF_phi;

  if (isum==0) {
    //printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET" rho(i,j,k)=%e\n",rho);
    //cout << "rho(i=0,j=0,k=0)=" << (*uu[0])[phi_index][theta_index][r_index][0] << endl;
   //printf("[init.cpp]: setup_averys_toyjet()\nrho(i,j,k)=%e\n",(*uu[0])[phi_index][theta_index][i][0]);
  //printf("[init.cpp]: setup_averys_toyjet()\n(i,j,k)=%d,%d,%d\n",i,j,k);
    //cout<<"[init.cpp]: rho[0][0][0]="<<*rho<<endl;
  }

  // int r_inspect=3; int theta_inspect=0; int phi_inspect=0;
  // if(r_index==r_inspect && theta_index==theta_inspect && phi_index==phi_inspect){
  //   printf(YELLOW"[setup_averys_toyjet.cpp]:"RED" (*uu[0])[%d][%d][%d][0]=%f\n"RESET,phi_index,theta_index,r_index,(*uu[0])[phi_index][theta_index][r_index][0]);
  //   printf(YELLOW"[setup_averys_toyjet.cpp]:"RED" (*uu[0])[%d][%d][%d][0]=%f, u_F_upt=%f,b_F_sq=%f,gamma=%f,gtt=%f,gtph=%f,OmegaOfPsi=%f\n"RESET,phi_index,theta_index,r_index,(*uu[0])[phi_index][theta_index][r_index][0],u_F_upt,b_F_sq,gamma,gtt,gtph,OmegaOfPsi);
  // }

  for (int var_index=0;var_index<=10;var_index++){
    if (var_index<0 || var_index >10) printf(YELLOW"[setup_averys_toyjet.cpp]: "RED"var_index=%d\n"RESET,var_index); // var_index had crazy values when I used (uu) instead of (*uu)

    //RG:FIXME
    //rest[var_index] = (*uu[0])[phi_index][theta_index][r_index][var_index];
    //rho=rest[0]*rhonor;              //density

    //RG:FIXME
    //rest[1]*=mp*cc*cc/3/kb/(rest[0]);//temperature
    // u[0]=rest[4];                    //4-velocity
    // u[1]=u[0]*rest[5];
    // u[2]=u[0]*rest[6];
    // u[3]=u[0]*rest[7];
    // Bi[1]=Bnor*rest[8];              //3-vector of magnetic field
    // Bi[2]=Bnor*rest[9];
    // Bi[3]=Bnor*rest[10];

    doub tmp=rest[var_index];

    if (isnan(tmp)) { //CHECK FOR NAN
      printf(YELLOW"[setup_averys_toyjet.cpp]: "RED"FOUND NAN IN ID (var_index=%d,r_idx=%d,th_idx=%d,phi_idx=%d)\n"RESET,var_index,r_idx,th_idx,phi_idx);
      
exit(1);
    }

  } //   for (int var_index=0;var_index<=10;var_index++){


  if (rest[0]<0) { //CHECK FOR UNPHYSICAL SOLUTIONS
    //printf(YELLOW"[setup_averys_toyjet.cpp]: "RED"FOUND UNPHYSICAL CONDITIONS IN ID (rest[0]=%f,r_idx=%d,th_idx=%d,phi_idx=%d)\n"RESET,rest[0],r_idx,th_idx,phi_idx);
    //exit(1);
  }

  //if ((th_idx==thlen/2) && (phi_idx==phlen/2)) printf(YELLOW"[setup_averys_toyjet.cpp]: "RESET"rho= %g r= %g\n",rest[0],r);

  // } //for r,theta,phi

  //if (isum==0) 
  //printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET" YO9\n");

  //exit(1);

  //return 0;
} //EOF





doub F_intp_lin(double F[2][rlen], double target_point)
{
  // INTERPOLATE F(Psi) to intp_point
  // y = y0 + (y1-y0)*(x-x0)/(x1-x0);
  
  // ASSUMES Psi positive
  
  //TODO: find nearest psi1,psi2 in F array 
  //target_point=5.2; // x
  // doub psi_lower=5.; // x0
  // doub psi_upper=5.4; // x1
  doub F_lower=0.; // y0
  doub F_upper=0.; // y1
  
  doub psi_lower = 0.; // just an initial value
  doub psi_upper = 1e10; // just an initial value
  doub dpsi=psi_upper-psi_lower; // x1-x0
  
  for (int r_idx=0; r_idx<rlen; r_idx++) { // LOOP THROUGH F(Psi) ARRAY
    if (F[0][r_idx] < target_point) {
      psi_lower = max(psi_lower,F[0][r_idx]); 
      F_lower   = F[1][r_idx];
    }
    if (F[0][r_idx] > target_point) {
      psi_upper = min(psi_upper,F[0][r_idx]); 
      F_upper   = F[1][r_idx];
    }
  }
  return ( F_lower + (F_upper-F_lower)*(target_point-psi_lower)/dpsi ); // y
}


void KS_metric(doub KS[4][4], doub iKS[4][4], doub r, doub costh)
{
  doub rsq=r*r;
  doub cossq=costh*costh;
  doub sinsq=1.-cossq;

  // RG: USED MATHEMATICA TO COMPUTE KERR-SCHILD METRIC AND ITS INVERSE
  
  // contravariant/inverse KS metric (see mathematica notebook)
  // CHECKED: THESE ARE INVERSE TO EACH OTHER TO MACHINE PRECISION
  // double KS[4][4];
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

  //-,+,+,+ signature covariant KS metric matrix (see e.g. eq(4) in McKinney & Gammie 2004)
  //RG: Why iKS? This is the *co*variant metric! iKS -> KS
  // double iKS[4][4];
  // iKS[0][0]=temp-1.;
  // iKS[0][1]=temp;
  // iKS[0][2]=0.;
  // iKS[0][3]=-a*temp*sinsq;
  // iKS[1][0]=iKS[0][1];
  // iKS[1][1]=1.+temp;
  // iKS[1][2]=0.;
  // iKS[1][3]=-a*(1.+temp)*sinsq;
  // iKS[2][0]=iKS[0][2];
  // iKS[2][1]=iKS[1][2];
  // iKS[2][2]=rhosq;
  // iKS[2][3]=0.;
  // iKS[3][0]=iKS[0][3];
  // iKS[3][1]=iKS[1][3];
  // iKS[3][2]=iKS[2][3];
  // iKS[3][3]=sinsq*(rhosq+asq*(1+temp)*sinsq);

}

void BL_metric(doub BL[4][4], doub iBL[4][4], doub r, doub costh)
{
  // Boyer-Lindquist metric
  // FIXME: BL.iBL != 1

  doub rsq=r*r;
  doub cossq=costh*costh;
  doub sinsq=1.-cossq;
  doub rhosq=rsq+asq*cossq; // akas Sigma
  doub Delta=rsq-2.*r+asq; // Lambda_squared in http://grwiki.physics.ncsu.edu/wiki/Kerr_Black_Hole
  doub temp=2.*r/rhosq;

  // if (isum==50) cout<<RED"ZERO OUT SPIN!"RESET<<endl; a=0.;asq=0.;

  // covariant (lower indices) BL metric matrix
  // see e.g. eq(3) in McKinney & Gammie 2004)
  // www.roma1.infn.it/teongrav/leonardo/bh/bhcap3.pdf
  // BL[0][0]=-(Delta/rhosq) - a*sinsq/rhosq; // does not look right...
  BL[0][0]=-(1. - 2.*r/rhosq); // eq (3.11)
  BL[0][1]=0.;
  BL[0][2]=0.;
  // BL[0][3]=Delta/rhosq * a *sinsq - sinsq/rhosq*(rsq+asq)*a; // does not look right...
  BL[0][3]=-2.*r/rhosq *a *sinsq; // eq (3.11)
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
  BL[3][3]=(rsq + asq + 2.*r*asq/rhosq*sinsq)*sinsq; // eq (3.11)

  // contravariant (upper indices) BL metric matrix
  // www.roma1.infn.it/teongrav/leonardo/bh/bhcap3.pdf
  iBL[0][0] = -(1./Delta) * ( rsq + asq + 2.*r*asq*sinsq/rhosq ); // eq (3.17)
  iBL[0][1] = 0.;
  iBL[0][2] = 0.;
  iBL[0][3] = -2.*r*a/Delta/rhosq; // eq (3.17)
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
}
