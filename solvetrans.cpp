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
  //RG:
  printf("[solvetrans] stNx=%d,sstep=%e\n",stNx,sstep);
	stNx++;                                 //increase count of points on each geodesic
	ittot++;                                //increase global count of points
    if(hh<sstep)                            //limit step size
		hh=sstep;
    int status = gsl_odeiv_evolve_apply (eee, ccc, ss, &sysx, &llog, ppy[currth].llmin, &hh, II);//evolve radiative transfer solution with a chosen step
    II[2]=fmod(II[2],(doub)2.*PI);          //phase winds up to be a large number, but only need phase modulo 2pi
    II[3]=fmod(II[3],(doub)2.*PI);
    if(II[0]<0)
		II[0]=0.;                           //plug for rare cases of negative intensity - such cases don't influence the final accuracy
    if(II[0]!=II[0]) {                      //check for NaNs
		printf("NaN in intensity evaluation, try smaller step size \n Exiting");
		exit(-1);
	}
};
gsl_odeiv_evolve_free (eee);                //free memory
gsl_odeiv_control_free (ccc);
gsl_odeiv_step_free (ss);
