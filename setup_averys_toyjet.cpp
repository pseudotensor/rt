//int setup_averys_toyjet(int i,int j,int k, float fluid_data[phlen][thlen][rlen][wdd]){ // i:r j:theta k:phi
{
  // Set input data for ASTRORAY (overwrite simulation data) 
  // to toy-jet+RIAF model by Broderick & Loeb 2009 
  // The Astrophysical Journal, 697:1164–1179, 2009 June 1
  // "IMAGING THE BLACK HOLE SILHOUETTE OF M87: 
  //  IMPLICATIONS FOR JET FORMATION AND BLACK HOLE SPIN"
  // http://dx.doi.org/10.1088/0004-637X/697/2/1164
  // http://adsabs.harvard.edu/abs/2009ApJ...697.1164B 
  // http://iopscience.iop.org/0004-637X/697/2/1164/ <- Reference for equation nrs
  // http://arxiv.org/abs/0812.0366
  // http://inspirehep.net/record/804072

// single fluid simulation dump files contains: 
// uu[?]:     [0    1   2     3                4    5    6        7      8    9        10   ]
//             rho, u, -u^t, -T^r_t/(rho u^r), u^t, v^r, v^theta, v^phi, B^r, B^theta, B^phi

  int r_index,theta_index,phi_index;
  // doub p_toyjet;
  for(theta_index=0;theta_index<thlen-1;theta_index++) // only boundaries
    for(phi_index=0;phi_index<phlen;phi_index++)
	  for(r_index=0;r_index<rlen;r_index++){
  
  int isum=r_index+theta_index+phi_index;
  
  if (isum==0) printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET" YO0\n");

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

  if (isum==0) printf("Get coordinates in here...\n");
  doub r     = coord[r_index][theta_index][0]; 
  doub theta = coord[r_index][theta_index][1];

  if (isum==0) cout <<"[init.cpp]: WIP What about dPhi here...?!"<<endl;
  doub phi   = ((doub) phi_index)/phlen / 2./PI; //RG: dPhi OFFSET?


  /*********************************
   * FREE PARAMETERS / USER CHOICE *
  /*********************************/
  doub xi_toyjet = 1.; // eq.(5) xi=0: cylindrical xi=1: conical determines collimation rate
  doub r_fp = 10.; // foot point of the jet where particles are injected
  doub p_toyjet = 2.-2.*xi_toyjet; // powerlaw-index in stream fct eq (A12) determines collimation rate 
  doub psi_toyjet = pow(r,p_toyjet) * (1. - cos(theta)); // stream fct eq.(5),(A42)
  doub OmegaOfPsi = pow(r,-3./2.);   // Kepler's law (floor@ISCO)
  //RG:FIXME:
  OmegaOfPsi = max(OmegaOfPsi , pow(6.,-3/2.));   // floor@ISCO


  // metric
  doub sigma=-1.; if (isum==0)printf("[CHECK:]sigma=%f\n",sigma); // sigma=u_F.u_F=+/-1 see between eq.(6),(7)

  if (isum==0) printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET"RG:FIXME HARDWIRED (flat) metric...\n");
  doub gtt=1.; doub gtphi=0.; doub gphiphi=r*r;
  doub detg=1.; if (isum==0)printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET"[HARDWIRE:]detg=%f\n",detg);
  if (isum==0) printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET" Need to adjust metric actually used in ASTRORAY...\n");


  // plasma
  // $ u^t_F $ eq.(A43)
  doub u_F_upt = 1./sqrt(gtt + 2.*gtphi + gphiphi*pow(OmegaOfPsi,2));
  doub u_F_upphi = u_F_upt*OmegaOfPsi; // follows from Omega definition just before eq.(A30)

  // $ b^\mu_F $ eq.(A44)
  doub b_F_upr     = -sigma * pow(r,p_toyjet)*sin(theta) / (u_F_upt * sqrt(-detg));
  doub b_F_uptheta = -sigma*p_toyjet * pow(r,p_toyjet-1)*(1.-cos(theta)) / (u_F_upt * sqrt(-detg));
  doub b_F_downphi = -2. * OmegaOfPsi * psi_toyjet * u_F_upt;
  doub b_F_downt   = -b_F_downphi * OmegaOfPsi;



  if (isum==0) printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET" FIXME: ...b_F_upt=?b_F_downt...\n"); doub b_F_upt=b_F_downt;
  if (isum==0) printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET" FIXME: ...b_F_upphi=?b_F_downphi...\n"); doub b_F_upphi=b_F_downphi;
  doub beta = sigma * b_F_upt / ( u_F_upt); // eq.(A46)
  if (isum==0) printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET" FIXME: ...b_F^2=?...\n");

  //RG:FIXME b_F: b_F^2 eq.(A40) ), u_F_upr, u_F_uptheta, u_F_uptheta 
  doub u_F_upr=0.; doub u_F_uptheta=0.; doub b_F_downr=0.; doub b_F_downtheta=0.; doub b_F_squared=0.;//RG:TODO

  // b_F_downr UNDEFINED, b_F_downtheta UNDEFINED, b_F_upphi UNDEFINED, U_F_upphi, b_F_upt
  // b_F_upphi I found out. see right before (A30)
  // b_F_upt comes from u_F^mu b_F_mu=0 see right after eq(8)
  doub b_F = sqrt(b_F_downr*b_F_upr + b_F_downtheta*b_F_uptheta + b_F_downphi/u_F_upt*(u_F_upt*b_F_upphi - u_F_upphi*b_F_upt)); // eq(A40)

  doub gamma = -sigma / sqrt( -(sigma + pow(beta*b_F,2) ) ); // eq.(A24)

  
  // JET (Sec. ?) & Appendix A

  // number density <-> free fct F  eq.(13)
  (*uu[0])[phi_index][theta_index][r_index][0]  = rhonor * u_F_upt * b_F_squared / gamma / (gtt + gtphi/OmegaOfPsi) * (1- exp(-pow(r/r_fp,2)/2.)); //eq.(13) F=1?

  // (*uu[0])[phi_index][theta_index][r_index][0]  = 1.;
  // if (isum==0) printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET" FIXME: ...density=1...\n");

  // n = rho_norm * exp(-pow(r/r_fp*sin(theta),2)/2.) // on that |r\cos\theta|=r_{fp} slice



  // float *u0     = (uu[0])[phi_index][theta_index][r_index][4]; // u^t? // u0 -> u_0 -> W
  // float *vr     = (uu[0])[phi_index][theta_index][r_index][5]; // Avery gives ur=u0*vr
  // float *vtheta = (uu[0])[phi_index][theta_index][r_index][6]; // Avery gives utheta=u0*vtheta
  // float *vphi   = (uu[0])[phi_index][theta_index][r_index][7]; // Avery gives uphi=u0*vphi

  // eq.(A46)
  (*uu[0])[phi_index][theta_index][r_index][4] = gamma*(u_F_upt + beta*b_F_upt);
  (*uu[0])[phi_index][theta_index][r_index][5] = gamma*(u_F_upr + beta*b_F_upr);
  (*uu[0])[phi_index][theta_index][r_index][6] = gamma*(u_F_uptheta + beta*b_F_uptheta);
  (*uu[0])[phi_index][theta_index][r_index][7] = gamma*(u_F_upphi + beta*b_F_upphi);


  // Magnetic field eq (A47)
  if (isum==0) printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET" FIXME: Compute b_F...\n");
  (*uu[0])[phi_index][theta_index][r_index][8]  = gamma * (b_F_upr - sigma*pow(b_F,2)*beta*u_F_upr);
  (*uu[0])[phi_index][theta_index][r_index][9]  = gamma * (b_F_uptheta - sigma*pow(b_F,2)*beta*u_F_uptheta);
  (*uu[0])[phi_index][theta_index][r_index][10] = gamma * (b_F_upphi - sigma*pow(b_F,2)*beta*u_F_upphi);

  // (*uu[0])[phi_index][theta_index][r_index][8]  = pow(r,p_toyjet-2.);
  // (*uu[0])[phi_index][theta_index][r_index][9]  = -p_toyjet * pow(r,p_toyjet-3.) * tan(theta/2.);
  // (*uu[0])[phi_index][theta_index][r_index][10] = -2*OmegaOfPsi * pow(r,p_toyjet-2.) * tan(theta/2.)/sin(theta);

  // SEG-faults
  // *Br     = (float)pow(r,p_toyjet-2.);
  // *Btheta = -p_toyjet * pow(r,p_toyjet-3.) * tan(theta/2.);
  // *Bphi   = -2*OmegaOfPsi * pow(r,p_toyjet-2.) * tan(theta/2.)/sin(theta);


  // DISK (RIAF) (Sec. ?)


  if (isum==0) {

    //RG: FIXME: PRINTS BUT THEN SEG_FAULTS...
    //cout << *rho << endl;
    printf(YELLOW"[init.cpp]:"RED" setup_averys_toyjet()\nrho(i,j,k)=%e\n"RESET,rho);

    //cout << "rho(i=0,j=0,k=0)=" << (*uu[0])[phi_index][theta_index][r_index][0] << endl;
   //printf("[init.cpp]: setup_averys_toyjet()\nrho(i,j,k)=%e\n",(*uu[0])[phi_index][theta_index][i][0]);
  //printf("[init.cpp]: setup_averys_toyjet()\n(i,j,k)=%d,%d,%d\n",i,j,k);


    //*rho = 42.;

    //cout<<"[init.cpp]: rho[0][0][0]="<<*rho<<endl;
  }


  //  return 0; // setup_averys_toyjet()

  if (isum==0) printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET" YO9\n");
      }
  }
