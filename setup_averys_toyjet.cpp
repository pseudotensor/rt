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

  //int r_index=ncut,theta_index=nt,phi_index=np;

  // doub p_toyjet;
  for(int theta_index=0;theta_index<thlen-1;theta_index++) // only boundaries
    for(int phi_index=0;phi_index<phlen;phi_index++)
      for(int r_index=0;r_index<rlen;r_index++){

        int isum=r_index+theta_index+phi_index;

        if (isum==0) printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET" Get coordinates in here...\n");
        //doub r     = exp(coord[r_index][theta_index][0]);
        doub r     = coord[r_index][theta_index][0]; 
        doub costh = coord[r_index][theta_index][1]; //RG:CHECK theta or cos(th)?
        rr=r;

        if (isum==0) cout <<YELLOW"[setup_averys_toyjet.cpp]:"RESET" WIP What about dPhi here...?!"<<endl;
        doub phi   = ((doub) phi_index)/phlen / 2./PI; //RG: dPhi OFFSET?
  
        costh=theta[r_index][theta_index];
        cossq=costh*costh;
        sinsq=1-cossq;
        sinth=sqrt(sinsq);
        rsq=rr*rr;
        rhosq=rsq+asq*cossq;
        Del=rsq-2.*rr+asq; // Lambda_squared in http://grwiki.physics.ncsu.edu/wiki/Kerr_Black_Hole                       
        doub temp=2.*rr/rhosq;

        // contravariant BL metric matrix                                                                                 
        iBL[0][0]=-(Del/rhosq) - a*sinsq/rhosq;
        iBL[0][1]=0.;
        iBL[0][2]=0.;
        iBL[0][3]=Del/rhosq * a *sinsq - sinsq/rhosq*(rsq+asq)*a;
        iBL[1][0]=iBL[0][1];
        iBL[1][1]=rhosq/Del;
        iBL[1][2]=0.;iBL[1][3]=0.;
        iBL[2][0]=iBL[0][2];
        iBL[2][1]=iBL[1][2];
        iBL[2][2]=rhosq;
        iBL[2][3]=0.;
        iBL[3][0]=iBL[0][3];
        iBL[3][1]=iBL[1][3];
        iBL[3][2]=iBL[2][3];
        iBL[3][3]=Del/rhosq*a*sinsq + sinsq/rhosq*(rsq+asq);


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



  /*********************************
   * FREE PARAMETERS / USER CHOICE *
  /*********************************/
  doub n0_th_RIAF  = 1.23e4; // 1.23x10^4 cm^-3 see main text after eq (4)
  doub n0_nth_RIAF = 3.8e2;  // 6.1 x10^2 cm^-3 see main text after eq (4)
  doub Te_RIAF=8.1e9;        // Te_0 = 8.1x10^9K in eq (2) see main text after eq (4)
  doub Rb=20; // Needed to reproduce for M87 (disk subdominant@7mm), see main text after eq (4)
  doub plasma_beta = 10;     // beta=10 see main text between eq (3) and (4)
  doub xi_toyjet = 1.; // eq.(5) xi=0: cylindrical xi=1: conical determines collimation rate
  doub r_fp = 10.; // foot point of the jet where particles are injected
  doub p_toyjet = 2.-2.*xi_toyjet; // powerlaw-index in stream fct eq (A12) determines collimation rate 
  doub psi_toyjet = pow(r,p_toyjet) * (1. - costh); // stream fct eq.(5),(A42)
  doub OmegaOfPsi = pow(r,-3./2.);   // Kepler's law (floor@ISCO)
  OmegaOfPsi = max(OmegaOfPsi , pow(6.,-3/2.));   // floor@ISCO

  // metric
  doub sigma=-1.; if (isum==0)printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET"[CHECK:]sigma=%f\n",sigma); // sigma=u_F.u_F=+/-1 see between eq.(6),(7)
  doub gtt=iBL[0][0];doub gtphi=iBL[0][3];doub gphiphi=iBL[1][1];
  doub grr=iBL[1][1];doub gthth=iBL[2][2];
  doub detg=-sinsq*pow(rsq+asq*cossq,2); // See e.g. arXiv:0706.0622v3 eq(70)

  // plasma
  // $ u^t_F $ eq.(A43)
  doub u_F_upt     = 1./sqrt( fabs(gtt + 2.*gtphi + gphiphi*pow(OmegaOfPsi,2)) );
  // JET:
  // BUT u_f.u_F=sigma -> u_F^t = Omega gtph/gtt +/- sqrt( (Omega gtph/gtt)^2 -Omega^2 gphph/gtt )
  doub u_F_upr     = 0.;                 //eq (6)
  doub u_F_uptheta = 0.;                 //eq (6)
  doub u_F_upphi   = u_F_upt*OmegaOfPsi; // follows from Omega definition just before eq.(A30)

  doub u_F_squared=sigma; //after eq.(A7)

  // $ b^\mu_F $ eq.(8,A44)
  doub b_F_downphi = -2. * OmegaOfPsi * psi_toyjet * u_F_upt;
  doub b_F_downt   = -b_F_downphi * OmegaOfPsi;

  // $ b^\mu_F $ eq.(7,A44) and some simple algebra
  doub b_F_upr     = -sigma * pow(r,p_toyjet)*sinth / (u_F_upt * sqrt(-detg));
  doub b_F_uptheta = -sigma*p_toyjet * pow(r,p_toyjet-1)*(1.-costh) / (u_F_upt * sqrt(-detg));
  doub b_F_upphi   = gtphi * b_F_downt + gphiphi * b_F_downphi;
  doub b_F_upt     = gtt * b_F_downt + gtphi * b_F_downphi;

  doub b_F_downr     = grr*b_F_upr;
  doub b_F_downtheta = gthth*b_F_uptheta;

  doub beta = sigma * b_F_upt / ( u_F_upt); // eq.(A46)
  // b_F^2 mostly cancels with beta^2 except for density...
  doub b_F_sq = b_F_downr*b_F_upr + b_F_downtheta*b_F_uptheta + b_F_downphi/u_F_upt*(u_F_upt*b_F_upphi - u_F_upphi*b_F_upt); // eq(A40)
  doub gamma = -sigma / sqrt( -(sigma + beta*beta*b_F_sq ) ); // eq.(A24)


  // u_F,b_F -> u,b : eqs (9,11,A46)
  doub u_RIAFtoyjet[4];
  u_RIAFtoyjet[0] = gamma * (u_F_upt     + beta*b_F_upt);
  u_RIAFtoyjet[1] = gamma * (u_F_upr     + beta*b_F_upr);
  u_RIAFtoyjet[2] = gamma * (u_F_uptheta + beta*b_F_uptheta);
  u_RIAFtoyjet[3] = gamma * (u_F_upphi   + beta*b_F_upphi);

  // DISK
  doub z=r*costh;
  doub Rc= sqrt(rsq-z*z); //cylindrical radius
  doub n_th_RIAF  = n0_th_RIAF  * exp(-z*z/(2.*Rc*Rc));
  doub n_nth_RIAF = n0_nth_RIAF * exp(-z*z/(2.*Rc*Rc));
  doub B_sq       = 8.*PI/plasma_beta *n_th_RIAF *mp*cc*cc/6./r; // eq(3) 
  doub B_RIAF_phi = sqrt(B_sq); // purely toroidal

  // number density RIAF (DISK)
  doub rho_RIAF = (n_th_RIAF + n_nth_RIAF)*mp;

  // JET (Sec. ?) & Appendix A <-> free fct F  eq.(13)
  rest[0]  = rho_RIAF + rhonor * u_F_upt * b_F_sq / gamma / (gtt + gtphi/OmegaOfPsi) * (1- exp(-pow(r/r_fp,2)/2.)); //eq.(13) F=1?
  rho/*[phi_index][theta_index][r_index]*/  = rho_RIAF + rhonor * u_F_upt * b_F_sq / gamma / (gtt + gtphi/OmegaOfPsi) * (1- exp(-pow(r/r_fp,2)/2.)); //eq.(13) F=1?
  (*uu[0])[phi_index][theta_index][r_index][0]  = rho_RIAF + rhonor * u_F_upt * b_F_sq / gamma / (gtt + gtphi/OmegaOfPsi) * (1- exp(-pow(r/r_fp,2)/2.)); //eq.(13) F=1?

  // (*uu[0])[phi_index][theta_index][r_index][0]  = 1.;
  // if (isum==0) printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET" FIXME: ...density=1...\n");

  // n = rho_norm * exp(-pow(r/r_fp*sin(theta),2)/2.) // on that |r\cos\theta|=r_{fp} slice



  // float *u0     = (uu[0])[phi_index][theta_index][r_index][4]; // u^t? // u0 -> u_0 -> W
  // float *vr     = (uu[0])[phi_index][theta_index][r_index][5]; // Avery gives ur=u0*vr
  // float *vtheta = (uu[0])[phi_index][theta_index][r_index][6]; // Avery gives utheta=u0*vtheta
  // float *vphi   = (uu[0])[phi_index][theta_index][r_index][7]; // Avery gives uphi=u0*vphi

  for(int u_idx=4;u_idx<=7;u_idx++) (*uu[0])[phi_index][theta_index][r_index][4] = u_RIAFtoyjet[u_idx];
  // (*uu[0])[phi_index][theta_index][r_index][4] = gamma*(u_F_upt + beta*b_F_upt);
  // (*uu[0])[phi_index][theta_index][r_index][5] = gamma*(u_F_upr + beta*b_F_upr);
  // (*uu[0])[phi_index][theta_index][r_index][6] = gamma*(u_F_uptheta + beta*b_F_uptheta);
  // (*uu[0])[phi_index][theta_index][r_index][7] = gamma*(u_F_upphi + beta*b_F_upphi);

  // Magnetic field eq (11,A47)
  (*uu[0])[phi_index][theta_index][r_index][8]  = gamma * (b_F_upr - sigma*b_F_sq*beta*u_F_upr);
  (*uu[0])[phi_index][theta_index][r_index][9]  = gamma * (b_F_uptheta - sigma*b_F_sq*beta*u_F_uptheta);
  (*uu[0])[phi_index][theta_index][r_index][10] = gamma * (b_F_upphi - sigma*b_F_sq*beta*u_F_upphi) + B_RIAF_phi;

  if (isum==0) {

    printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET" rho(i,j,k)=%e\n",rho);

    //cout << "rho(i=0,j=0,k=0)=" << (*uu[0])[phi_index][theta_index][r_index][0] << endl;
   //printf("[init.cpp]: setup_averys_toyjet()\nrho(i,j,k)=%e\n",(*uu[0])[phi_index][theta_index][i][0]);
  //printf("[init.cpp]: setup_averys_toyjet()\n(i,j,k)=%d,%d,%d\n",i,j,k);
    //cout<<"[init.cpp]: rho[0][0][0]="<<*rho<<endl;
  }

  //  return 0; // setup_averys_toyjet()

  if (isum==0) printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET" (*uu[0])[0][0][0][0]=%f\n",(*uu[0])[phi_index][theta_index][r_index][0]);


  // for(m=0;m<11;m++) rest[m]=(*uu[fdiff])[np][nt][xrn][m];

  for (int var_index=0;var_index<=10;var_index++){
    if (var_index<0 || var_index >10) printf(YELLOW"[setup_averys_toyjet.cpp]: "RED"var_index=%d\n"RESET,var_index);
    rest[var_index] = (*uu[0])[phi_index][theta_index][r_index][var_index];
    rho=rest[0]*rhonor;              //density
    rest[1]*=mp*cc*cc/3/kb/(rest[0]);//temperature
    u[0]=rest[4];                    //4-velocity
    u[1]=u[0]*rest[5];
    u[2]=u[0]*rest[6];
    u[3]=u[0]*rest[7];
    Bi[1]=Bnor*rest[8];              //3-vector of magnetic field
    Bi[2]=Bnor*rest[9];
    Bi[3]=Bnor*rest[10];

    float tmp=(*uu[0])[phi_index][theta_index][r_index][var_index];
    //if (isum==0) printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET" tmp=%f\n",tmp);

    if (phi_index==0 && theta_index==0 && r_index==10) { // nan in rho,8 & 10
      printf(YELLOW"[setup_averys_toyjet.cpp]: "RED"tmp=%f var_index=%d\n"RESET,tmp,var_index);
      printf(YELLOW"[setup_averys_toyjet.cpp]: "RED"(*uu[0])[phi_index][theta_index][r_index][var_index]=%f var_index=%d\n"RESET,(*uu[0])[phi_index][theta_index][r_index][var_index],var_index);
    }

    // if (tmp!=tmp) { //CHECK FOR NAN // unsafe in -Ofast or -ffast-math
    if (isnan(tmp)) { //CHECK FOR NAN
      printf(YELLOW"[setup_averys_toyjet.cpp]: "RED"FOUND NAN IN ID (var_index=%d,r_index=%d,th_index=%d,ph_index=%d)\n"RESET,var_index,r_index,theta_index,phi_index);
      exit(1);
    }

  } //   for (int var_index=0;var_index<=10;var_index++){

  // exit(1);

  } //for r,theta,phi

  printf(YELLOW"[setup_averys_toyjet.cpp]: "RED"(*uu[0])[0][0][10][0]=%f\n"RESET,(*uu[0])[0][0][10][0]); // -nan why does nan checker above does not catch it?!

  //if (isum==0) 
  printf(YELLOW"[setup_averys_toyjet.cpp]:"RESET" YO9\n");

} //EOF
