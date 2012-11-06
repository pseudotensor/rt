//void solvetrans (void)
doub II[5], hh,llog;int stNx;//int k; //abs accuracy 1e-6 and method 
gsl_odeiv_step *ss; gsl_odeiv_control *ccc; gsl_odeiv_evolve *eee;
gsl_odeiv_system sysx = {trans, NULL, 5, &currth};
eee = gsl_odeiv_evolve_alloc(5);
ccc = gsl_odeiv_control_standard_new(accurr, accurr, 1.0, 0.0);
ss = gsl_odeiv_step_alloc (TT, 5);

for(k=0;k<5;k++)II[k]=Iint;stNx=0;
llog=ppy[currth].llmax;
II[2]=Iang;II[3]=II[2];
hh=sstep;
//if(llog<5){int hj=1;}

while (llog > ppy[currth].llmin-1.3*sstep)
{	stNx++;ittot++;
if(hh<sstep)hh=sstep;
int status = gsl_odeiv_evolve_apply (eee, ccc, ss, &sysx, &llog, ppy[currth].llmin, &hh, II);
II[2]=fmod(II[2],(doub)2.*PI);
II[3]=fmod(II[3],(doub)2.*PI);

//int status = gsl_odeiv_step_apply (ss, llog, hh,II, yerr, NULL, NULL, &sysx);llog+=hh;

if(II[0]<0)II[0]=0.;//if(II[4]<0)II[4]=0.;
if(II[0]>0.1)II[0]=0.;//cuts run-away errors
if(II[0]!=II[0]){hh=sstep;};
//if(II[4]>100000.*II[0]){hh=step;}
;};//printf("Done %f\n",llog);
//printf("N=%d, stNx=%d\n",omp_get_thread_num(),stNx);
gsl_odeiv_evolve_free (eee);gsl_odeiv_control_free (ccc);gsl_odeiv_step_free (ss);