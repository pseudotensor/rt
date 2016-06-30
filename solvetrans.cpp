time_t t_b4_solvetrans = clock();

const int nr_of_imagediags=5;
doub II[nr_of_imagediags],                                 //array of intensities
	 hh,                                    //step of radiative transfer solver
	 llog;                                  //log-parameter on a geodesic
int  k,
	 counter_pt_on_geodesic=0;                                //counter of radiative transfer point along geodesic
gsl_odeiv_step *gsl_aux_ode_vars;                         //auxiliary ODE variables from GSL (GNU Scientific Library)
gsl_odeiv_control *gsl_accuracy; 
gsl_odeiv_evolve *gsl_alloc;
const gsl_odeiv_step_type * gsl_time_evol;             //ODE solver type
gsl_time_evol = gsl_odeiv_step_rk2;                    //use Runge-Kutte 2-nd order

//DEFAULT:
gsl_odeiv_system trans_gsl = {trans, NULL, nr_of_imagediags, &currth};//prepare to solve ODE of nr_of_imagediags variables given by "trans" function
//RG:WIP
// gsl_odeiv_system trans_gsl = {trans, NULL, nr_of_imagediags, &currth, counter_pt_on_geodesic};//prepare to solve ODE of nr_of_imagediags variables given by "trans" function

gsl_alloc = gsl_odeiv_evolve_alloc(nr_of_imagediags);
gsl_accuracy = gsl_odeiv_control_standard_new(accurr, accurr, 1.0, 0.0);//relative AND absolute accuracy accurr
gsl_aux_ode_vars = gsl_odeiv_step_alloc (gsl_time_evol, nr_of_imagediags);

for(k=0;k<nr_of_imagediags;k++)                           //initialize intensities at small values (can't initialize with zeros)
	II[k]=Iint;
II[2]=Iang;                                //since variables 2 and 3 represent angles in "polarized coordinates", they are given values order 1
II[3]=II[2];
llog=ppy[currth].llmax;                    //maximum log-parameter
hh=sstep;                                  //initial step (step is always negative)

while (llog > ppy[currth].llmin-1.3*sstep) { // MAIN RADIATIVE TRANSFER LOOP
	counter_pt_on_geodesic++;                // increase count of points on each geodesic
	ittot++;                                 // increase global count of points

    // original ASTRORAYv1.0 BUT: RG:FIXME sstep<0 -> floor step by hh>sstep
    if(hh<sstep) hh=sstep;                  //floor on step size
    // if(hh>sstep) hh=sstep;                  //floor on step size //GOT LOOKUP ERROR OF CLOSEST GEODESIC POINT [evalpointzero.cpp]
    //if (counter_pt_on_geodesic==1) printf(YELLOW"[solvetrans.cpp]: "RESET"Calling trans()... sstep=hh=%f\nCHANGED hh<sstep -> hh>sstep\n2DO:CATCH CASE WHEN GETTING STUCK IN trans()\n",hh);

    // actual call to trans() fct: solve radiative transfer equation
    //                                   vars,acc,    trans()
    int status = gsl_odeiv_evolve_apply (gsl_alloc, gsl_accuracy, gsl_aux_ode_vars, &trans_gsl, 
    // &llog, ppy[currth].llmin, &hh, II, counter_pt_on_geodesic); //FCT ARG TO TRANS() //RG:WIP
                                         &llog, ppy[currth].llmin, &hh, II); // FCT ARG TO TRANS()
    // doub tet = gsl_odeiv_evolve_apply (gsl_alloc, gsl_accuracy, gsl_aux_ode_vars, &trans_gsl, 
    //                                      // &llog, ppy[currth].llmin, &hh, II, counter_pt_on_geodesic); //FCT ARG TO TRANS() //RG:WIP
    //                                      &llog, ppy[currth].llmin, &hh, II); // FCT ARG TO TRANS()


    II[2]=fmod(II[2],(doub)2.*PI);          //phase winds up to be a large number, but only need phase modulo 2pi
    II[3]=fmod(II[3],(doub)2.*PI);

    if (II[0]<0) { // I<0
      printf(YELLOW"[solvetrans.cpp]: "RED"Negative intensity, maybe try smaller step size!?\nWill continue...\n"RESET);
	  II[0]=0.;                           //plug for rare cases of negative intensity - such cases don't influence the final accuracy
    }

    if(isnan(II[0])) {                      //check for NaNs
		printf("NaN in intensity evaluation, try smaller step size \n Exiting");
        t_solvetrans += (clock() - t_b4_solvetrans) / (doub)CLOCKS_PER_SEC;
		exit(-1);
	}

    //RG:FIXME WORKS HERE BUT SLOWS THE CODE DOWN ENORMOUSLY
    //t_solvetrans += (clock() - t_b4_solvetrans) / (doub)CLOCKS_PER_SEC;
    if ( (counter_pt_on_geodesic==1) && (currth==0) ) {
      // if (false) {   
    if (TEMPERATURE_DIAGNOSTIC) {   
      /****************
       * TEMPERATURES *
       * given T_sim  *
       * get   Te,Tp  *
       ****************/

      stringstream snapshot;
       
      snapshot<<(int)fnum;
      ofstream temperature_harm_grid_file;
      temperature_harm_grid_file.open((dir+astr[sp]+xstr+"T_grid_"+astr[sp]+"-snapshot"+snapshot.str()+".dat").c_str(), ios::out | ios::app | ios::binary);

      doub rest[11];

  for(int r_idx=0;r_idx<ndd;r_idx++) // rlen ~> rcut ?
    for(int th_idx=0;th_idx<thlen;th_idx++)
      for(int ph_idx=0;ph_idx<phlen;ph_idx++) {

        doub rr=exp(coord[r_idx][th_idx][0]); // rtab[r_idx]; // rlen ~> rcut ? Likely wrong for r>rcut !
        doub costh=coord[r_idx][th_idx][1];

      //RG:TODO get coordinates
      doub ph;
      // get_coordinates(rr, costh, ph, currth, counter_pt_on_geodesic);
      doub tet;
      doub tpt,T_sim,rho; 
      //RG:TODO get T_sim, rho
      // get_fluid(rr, costh, ph, T_sim, rho);
      
      if (rr <= rcut) {  // INSIDE SIMULATION ZONE
        for(int m=0;m<11;m++) rest[m]=(*uu[0])[ph_idx][th_idx][r_idx][m];

        T_sim=rest[1]*mp*cc*cc/3/kb/rest[0];//internal energy density

      } else { // OUTSIDE SIMULATION ZONE: RADIAL EXTENSION
        // for(int m=0;m<5;m++) rest[m]=uext[ph_idx][th_idx][m]; // using quantites at sphere r=rcut for radii>rcut for power-law extension

        for(int m=0;m<11;m++) rest[m]=(*uu[fdiff])[ph_idx][th_idx][ncut-1][m];
		rho=rest[0]*rhonor;               // density
		rest[1]*=mp*cc*cc/3/kb/(rest[0]); // temperature //RG:CHECK UNITS WHEN USING AVERY'S MODEL
        uext[ph_idx][th_idx][4]=rest[1];

        // T_sim=uext[ph_idx][th_idx][4]*pow(rr/rcut,-1.0);   //temperature is extended with power-law index "-1.0"
        T_sim=uext[ph_idx][th_idx][4]*pow(rr/rcut,-Upo);   //temperature is extended with power-law index "-1.0"

      } // if (rr < rcut) {  // INSIDE SIMULATION ZONE

      get_electron_temperature (T_sim, 0., tet, tpt);
      // temperature_diag(rr, costh, ph, T_sim, tet, tpt, rho, rhonor, currth, counter_pt_on_geodesic);

      if (r_idx+th_idx+ph_idx==0) 
        // temperature_harm_grid_file<<"# r_idx th_idx ph_idx T_e T_p T_sim rho u r costh"<<endl; // WRITE HEADER
        temperature_harm_grid_file<<"# r_idx th_idx ph_idx T_e T_p T_sim rho r costh"<<endl; // WRITE HEADER

      temperature_harm_grid_file<<r_idx<<" "<<th_idx<<" "<<ph_idx<<" "<<tet<<" "<<tpt<<" "<<T_sim<<" "<<rho/*<<" "<<u<<" "<<magn*/<<" "<<rr<<" "<<costh<<endl; // can add tet_isoth etc later        
 
      } // for (r,th,ph)
  temperature_harm_grid_file.close();

    } // if (TEMPERATURE_DIAGNOSTIC) {   
    } // if ( (counter_pt_on_geodesic==0) && (currth==0) ) {

    // exit(1);

    // TODO: only one geodesic
    // if ( (counter_pt_on_geodesic==0) && (rr>=0.99*exp(coord[ncuttab[0]-1][0][0])) && (rr<=1.01*exp(coord[ncuttab[0]-1][0][0])) || (rr>=0.99*exp(coord[ncuttab[0]][0][0])) && (rr<=1.01*exp(coord[ncuttab[0]][0][0])) || (rr>=0.99*exp(coord[ncuttab[0]+1][0][0])) && (rr<=1.01*exp(coord[ncuttab[0]+1][0][0])) ) {
    //   printf(YELLOW"[evalpointzero.cpp]: "MAGENTA"T_sim=%g,tet=%g,tpt=%g,rr=%g"RESET"\n",T_sim,tet,tpt,rr);
    //   // exit(-1);
    // }


 }; // while (llog > ppy[currth].llmin-1.3*sstep) { // MAIN RADIATIVE TRANSFER LOOP

gsl_odeiv_evolve_free (gsl_alloc);                //free memory
gsl_odeiv_control_free (gsl_accuracy);
gsl_odeiv_step_free (gsl_aux_ode_vars);

//RG:FIXME t_solvetrans=0 :-S
//t_solvetrans += (clock() - t_b4_solvetrans) / (doub)CLOCKS_PER_SEC;
