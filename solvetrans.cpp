time_t t_b4_solvetrans = clock();

doub II[5],                                 //array of intensities
	 hh,                                    //step of radiative transfer solver
	 llog;                                  //log-parameter on a geodesic
int  k,
	 stNx=0;                                //counter of radiative transfer point along geodesic
gsl_odeiv_step *ss;                         //auxiliary ODE variables from GSL (GNU Scientific Library)
gsl_odeiv_control *ccc; 
gsl_odeiv_evolve *eee;
const gsl_odeiv_step_type * TT;             //ODE solver type
TT = gsl_odeiv_step_rk2;                    //use Runge-Kutte 2-nd order
gsl_odeiv_system sysx = {trans, NULL, 5, &currth};//prepare to solve ODE of 5 variables given by "trans" function
eee = gsl_odeiv_evolve_alloc(5);
ccc = gsl_odeiv_control_standard_new(accurr, accurr, 1.0, 0.0);//relative AND absolute accuracy accurr
ss = gsl_odeiv_step_alloc (TT, 5);

for(k=0;k<5;k++)                           //initialize intensities at small values (can't initialize with zeros)
	II[k]=Iint;
II[2]=Iang;                                //since variables 2 and 3 represent angles in "polarized coordinates", they are given values order 1
II[3]=II[2];
llog=ppy[currth].llmax;                    //maximum log-parameter
hh=sstep;                                  //initial step (step is always negative)

while (llog > ppy[currth].llmin-1.3*sstep) {//main radiative transfer loop
	stNx++;                                 //increase count of points on each geodesic
	ittot++;                                //increase global count of points

    //RG:FIXME sstep<0 -> floor step by hh>sstep
    if(hh<sstep) hh=sstep;                  //floor on step size
    // if(hh>sstep) hh=sstep;                  //floor on step size //GOT LOOKUP ERROR OF CLOSEST GEODESIC POINT [evalpointzero.cpp]

    //if (stNx==1) printf(YELLOW"[solvetrans.cpp]: "RESET"Calling trans()... sstep=hh=%f\nCHANGED hh<sstep -> hh>sstep\n2DO:CATCH CASE WHEN GETTING STUCK IN trans()\n",hh);

    // actual call to trans() fct: solve radiative transfer equation
    //                                   vars,acc,    trans()
    int status = gsl_odeiv_evolve_apply (eee, ccc, ss, &sysx, 
                                         // &llog, ppy[currth].llmin, &hh, II, stNx); //FCT ARG TO TRANS() //RG:WIP
                                         &llog, ppy[currth].llmin, &hh, II); //FCT ARG TO TRANS()

    II[2]=fmod(II[2],(doub)2.*PI);          //phase winds up to be a large number, but only need phase modulo 2pi
    II[3]=fmod(II[3],(doub)2.*PI);

    if (II[0]<0) { // I<0
      //printf(YELLOW"[solvetrans.cpp]: "RED"Negative intensity, maybe try smaller step size!?\nWill continue...\n"RESET);
	  II[0]=0.;                           //plug for rare cases of negative intensity - such cases don't influence the final accuracy
    }

    if(isnan(II[0])) {                      //check for NaNs
		printf("NaN in intensity evaluation, try smaller step size \n Exiting");
        t_solvetrans += (clock() - t_b4_solvetrans) / (doub)CLOCKS_PER_SEC;
		exit(-1);
	}

    //RG:FIXME WORKS HERE BUT SLOWS THE CODE DOWN ENORMOUSLY
    //t_solvetrans += (clock() - t_b4_solvetrans) / (doub)CLOCKS_PER_SEC;
};

gsl_odeiv_evolve_free (eee);                //free memory
gsl_odeiv_control_free (ccc);
gsl_odeiv_step_free (ss);

//RG:FIXME t_solvetrans=0 :-S
//t_solvetrans += (clock() - t_b4_solvetrans) / (doub)CLOCKS_PER_SEC;
