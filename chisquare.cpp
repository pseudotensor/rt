// void chisquare(doub F[], doub LP[], doub CP[], doub chisq, doub chisq_I, int dof=12, int dof_I=7, int n_par=3)
void chisquare(doub F[], doub LP[], doub CP[], doub& chisq, doub& chisq_I, doub (&resid)[14],int dof=12, int dof_I=7, int n_par=3, string CHOOSE_CHI=DRIVE_FIT_TO) //dof=12 ~> include polarization
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

  // doub resid[sflen];

  // [nu, F, LP, EVPA, CP]
  const doub tofit[sflen][5]={
  {8.450, 0.683, 0., 0., -0.2500},
  {14.90, 0.871, 0., 0., -0.6200},
  {22.50, 0.979, 0.1900, 131.0, 0.},
  {43.00, 1.135, 0.5500, 94.25, 0.},
  {87.73, 1.841, 1.420, -4., 0.},
  {102.0, 1.908, 0., 0., 0.},
  {145.0, 2.275, 0., 0., 0.},
  //  {230.9, 2.637, 7.398, 111.5, -1.200}, // ASTRORAYv1.0
  {230.9, 2.637, 6., 111.5, -1.200}, // LP probably more like 5% ...
  {349.0, 3.181, 6.499, 146.9, -1.500},
  {674.0, 3.286, 0., 0., 0.},
  {857.0, 2.867, 0., 0., 0.},
  {1500., 1., 0., 0., 0.},{3000., 1., 0., 0., 0.},{5000., 1., 0., 0., 0.}};

  doub F230GHz_Dexter=3.75, dF230Ghz_Dexter=0.48, spectral_index=-0.18, dspectral_index=0.34; // Dexter+2010 Submillimeter bump paper
  doub chisq_Dexter=pow((F[7]-F230GHz_Dexter)/dF230Ghz_Dexter,2) + pow(((F[9]-F[7])/(tofit[9][0]-tofit[7][0]) - spectral_index)/dspectral_index,2); // range for slope is more like [230GHz-770GHz]

  // CK's chi^2 ...WIP... chi^2=((sim - obs) / sigma)^2
  // nu    = [95.589d, 151.77d, 359.224d, 660.735d] * 1d9
  // obs   = [ 1.897d, 2.611d, 3.592d, 2.688d]
  // sigma = [0.639+0.645,0.878+0.872,1.303+1.267,1.385+1.335] / 2d
  doub F_CK[4][3] = { // {nu, F, dF}
    {95.589, 1.897, 0.64}, 
    {151.77, 2.611, 0.9}, 
    {359.224, 3.592, 1.3},
    {660.735, 2.688, 1.3}
  };
  doub chisq_CK=0;
  // for( int CK_idx=0; CK_idx<=3; CK_idx++) {
  //   chisq_CK+= (F[]-F_CK[CK_idx][1])/F_CK[CK_idx][2]; // interpolate or change modelled frequencies ?
  // }
  chisq_CK+= pow((F[5]-F_CK[0][1])/F_CK[0][2],2); // frequencies inconsistent nu=102 vs nu=96
  chisq_CK+= pow((F[6]-F_CK[1][1])/F_CK[1][2],2); // frequencies inconsistent nu=145 vs nu=152
  chisq_CK+= pow((F[8]-F_CK[2][1])/F_CK[2][2],2); // frequencies inconsistent nu=359 vs nu=349
  chisq_CK+= pow((F[9]-F_CK[3][1])/F_CK[3][2],2); // frequencies inconsistent nu=674 vs nu=661

  printf(YELLOW"[chisquare.cpp]: "RESET"chi^2_Dexter=%f chi^2_CK=%f\n",chisq_Dexter,chisq_CK);

  // measurement errors of mean fluxes, CP fractions, LP fractions, and EVPAs (no measurements at highest frequencies)                                     
  const doub dFnu[sflen]={0.031, 0.012, 0.015, 0.026, 0.080, 0.1517, 0.2644, 0.1414, 0.1205, 0.3508, 0.2404, 0., 0., 0.},
    dCP=0.30,                    // at 230GHz and 345GHz              
      // dLP[3]={0.50, 0.658, 0.605}, // at 87GHz, 230GHz, and 345GHz
    dLP[3]={1., 1., 1.}, // at 87GHz, 230GHz, and 345GHz
    dEVPA[3]={11.,5.4,2.21};     // at 87GHz, 230GHz, and 345GHz

  const bool trustLP87=true;//whether to fit for LP fraction at 87GHz. Its observed value is controversial
 

      /***********************************************************************************/
      if (DRIVE_FIT_TO=="ASTRORAYv1.0") {

        for(int il=4;il<=10;il++)           // LOOP THROUGH nu only up until k<=10 because thats where there are errorbars // and check for kmin>=4
          resid[il-4]=(F[il]-tofit[il][1])/dFnu[il]; // F         

        if (trustLP87) resid[7]=(LP[4]-tofit[4][2])/dLP[0]; // LP nu= 87Ghz
        else resid[7] = 0.;

        // resid[7]=(LP[4]-tofit[4][2])/dLP[0]; // LP nu= 87Ghz
        resid[8]=(LP[7]-tofit[7][2])/dLP[1]; // LP nu=230Ghz
        resid[9]=(LP[8]-tofit[8][2])/dLP[2]; // LP nu=345Ghz
        resid[10]=(CP[7]-tofit[7][4])/dCP;    // CP nu=230Ghz
        resid[11]=(CP[8]-tofit[8][4])/dCP;    // CP nu=345Ghz

      }
      else if (DRIVE_FIT_TO=="ASTRORAYv1.0_unpolarized") {

        for(int il=4;il<=10;il++)           // LOOP THROUGH nu only up until k<=10 because thats where there are errorbars // and check for kmin>=4
          resid[il-4]=(F[il]-tofit[il][1])/dFnu[il]; // F 
      }
      else if (DRIVE_FIT_TO=="DEXTER") {
        
      }
      else {
        printf(YELLOW"[m_sear.cpp]: "RESET"YOU NEED TO SPECIFY THE \"DRIVE_FIT_TO\" VARIABLE in [global_variables.cpp]\n");
        exit(1);
      }
      // RG: COULD ADD MORE (New Flux measurements,EVPA,source size,...) HERE ...WIP...
      /***********************************************************************************/

   
  if (dof>7) { // WANT2 INCLUDE POLARIZATION?
  }

  for(int il=0;il<dof;il++) {

    // WIP: REARRANGE resid... incorporate I,LP,CP differently...
    chisq+=resid[il]*resid[il]/(dof-n_par); //calculate \chi^2/dof , 3 parameters are varied: inclination, rho_nor, C
    if (il<dof_I) chisq_I+=resid[il]*resid[il]/(dof_I-n_par); //calculate \chi^2/dof , 3 parameters are varied: inclination, rho_nor, C
  //RG:FLAG SHOULD INCORPORATE
  //if(trustLP87)            //either include or not include in the fit the LP fraction at 87GHz
  }

  // talk about it
  // printf("chisquare(): chi^2/dof=%f chi^2_I/dof=%f\n",chisq,chisq_I);
  for(int kk=4;kk<=10;kk++){
    printf(YELLOW"[chisquare.cpp]:"RESET"  avg at f=%.1f; I=%.3fJy I_obs=%.3fJy LP=%.2f%% LP_obs=%.2f%% CP=%.3f%% CP_obs=%.3f%% EVPA=%.2f \n",sftab[kk][0], F[kk],tofit[kk][1],LP[kk],tofit[kk][2],CP[kk],tofit[kk][4],xEVPA[kk]);
  }
  
};
