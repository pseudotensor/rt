// void chisquare(doub F[], doub LP[], doub CP[], doub chisq, doub chisq_I, int dof=12, int dof_I=7, int n_par=3)
void chisquare(doub F[], doub LP[], doub CP[], doub& chisq, doub& chisq_I, int dof=12, int dof_I=7, int n_par=3) //dof=12 ~> include polarization
{
  /***********************************************************************
   * GIVEN FLUXES F, LINEAR POLARIZATION LP AND CIRCULAR POLARIZATION CP *
   * COMPUTE CHI^2 (polarized: "chisq" and unpolarized "chisq_I")        *
   ***********************************************************************/

  const int sflen=14; // NR OF MEASURED FLUXES
  // int dof_I= 7; // degrees of freedom (total I only)
  // int dof  =12; // degrees of freedom (I+LP+CP)
  // int n_par= 3; // nr of free parameters we vary (th,heat,rhonor)

  // CONSTRAINTS FROM OBSERVATIONS, SEE SECTION 2 AND TABLE I IN:
  // http://adsabs.harvard.edu/abs/2012ApJ...755..133S)
  // http://arxiv.org/abs/1007.4832

  doub resid[sflen];

  // [nu, F, LP, EVPA, CP]
  const doub tofit[sflen][5]={
  {8.450, 0.683, 0., 0., -0.2500},
  {14.90, 0.871, 0., 0., -0.6200},
  {22.50, 0.979, 0.1900, 131.0, 0.},
  {43.00, 1.135, 0.5500, 94.25, 0.},
  {87.73, 1.841, 1.420, -4., 0.},
  {102.0, 1.908, 0., 0., 0.},
  {145.0, 2.275, 0., 0., 0.},
  {230.9, 2.637, 7.398, 111.5, -1.200}, // LP probably more like 5% ...
  {349.0, 3.181, 6.499, 146.9, -1.500},
  {674.0, 3.286, 0., 0., 0.},
  {857.0, 2.867, 0., 0., 0.},
  {1500., 1., 0., 0., 0.},{3000., 1., 0., 0., 0.},{5000., 1., 0., 0., 0.}};

  // measurement errors of mean fluxes, CP fractions, LP fractions, and EVPAs (no measurements at highest frequencies)                                     
  const doub dFnu[sflen]={0.031, 0.012, 0.015, 0.026, 0.080, 0.1517, 0.2644, 0.1414, 0.1205, 0.3508, 0.2404, 0., 0., 0.},
    dCP=0.30,                    // at 230GHz and 345GHz              
    dLP[3]={0.50, 0.658, 0.605}, // at 87GHz, 230GHz, and 345GHz      
    dEVPA[3]={11.,5.4,2.21};     // at 87GHz, 230GHz, and 345GHz      

  const bool trustLP87=true;//whether to fit for LP fraction at 87GHz. Its observed value i s controversial
    
  for(int il=4;il<=10;il++) {          // LOOP THROUGH nu only up until k<=10 because thats where there are errorbars // and check for kmin>=4
	resid[il-4]=(F[il]-tofit[il][1])/dFnu[il]; // F         
  }
  if (dof>7) { // WANT2 INCLUDE POLARIZATION?
    resid[7]=(LP[4]-tofit[4][2])/dLP[0]; // LP nu= 87Ghz
    resid[8]=(LP[7]-tofit[7][2])/dLP[1]; // LP nu=230Ghz
    resid[9]=(LP[8]-tofit[8][2])/dLP[2]; // LP nu=345Ghz
    resid[10]=(CP[7]-tofit[7][4])/dCP;    // CP nu=230Ghz
    resid[11]=(CP[8]-tofit[8][4])/dCP;    // CP nu=345Ghz
  }


  // doub chisq=0.,chisq_I=0.;

  for(int il=0;il<dof;il++) {

    // WIP: REARRANGE resid... incorporate I,LP,CP differently...
    chisq+=resid[il]*resid[il]/(dof-n_par); //calculate \chi^2/dof , 3 parameters are varied: inclination, rho_nor, C
    if (il<dof_I) chisq_I+=resid[il]*resid[il]/(dof_I-n_par); //calculate \chi^2/dof , 3 parameters are varied: inclination, rho_nor, C
  //RG:FLAG SHOULD INCORPORATE
  //if(trustLP87)            //either include or not include in the fit the LP fraction at 87GHz
  }

  // talk about it
  // printf("chisquare(): chi^2/dof=%f chi^2_I/dof=%f\n",chisq,chisq_I);
  
};
