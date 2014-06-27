int p,                     //for loop variable
	stN=0;                 //counter of points along geodesic
doub tpast,                //time at a previous step
	 y[12],                //main array for a geodesic solver
	 der[12],              //auxiliary array of derivates, invoked by ODE solver
	 txx,                  //auxiliary affine parameter for sinh conversion (see below)
	 lam[maxco],           //affine parameter along geodesic
	 coox[12][maxco];      //coordinates, tangential vector, and perpendicular parallel-tranpsported vector in each point along geodesic
gsl_odeiv_step *s;         //auxiliary ODE variables from GSL (GNU Scientific Library)
gsl_odeiv_control *c; 
gsl_odeiv_evolve *e;

int currth = omp_get_thread_num(); //OpenMP thread ID => geodesic number
const gsl_odeiv_step_type * T;     //ODE solver type
T = gsl_odeiv_step_rk2;            //use Runge-Kutte 2-nd order
s = gsl_odeiv_step_alloc (T, 12);  //initialize with 12 variables
c = gsl_odeiv_control_standard_new(0.0, accur, 1.0, 0.0);//relative accuracy "accur"
e = gsl_odeiv_evolve_alloc (12);
gsl_odeiv_system sys = {solvegeodesic, NULL, 12, NULL};//prepare to solve ODE of 12 variables given by "solvegeodesic" function

//a - spin, r0 - distance from picture plain to BH, th - cosine of polar angle, b - impact parameter, 
//beta - angle in a picture plane counterclockwise from the direction to the north pole  
//BH rotates clockwise as viewed from above the north pole

rg=1+sqrt(1.-asq);             //BH horizon
doub bsq=b*b,                  //auxiliary geometric variables based on a, r0, th, b
	 sinbeta=sin(beta),
	 cosbeta=cos(beta),
	 sinbetasq=sinbeta*sinbeta,
	 cosbetasq=1-sinbetasq,
	 costh=cos(th),
	 sinthsq=1-costh*costh,
	 sinth=sqrt(sinthsq), 
	 costhsq=1-sinthsq, 
	 r0sq=r0*r0,
	 rmin2=r0-2.;
//initial quantities in the picture plane
doub tin = 0.,                                                                            //initial time t(0)
	 rin = 0.5*bsq/r0 + r0,                                                               //r(0) - not precisel r0 because of curvature, but assumed b<<r0
	 phiin = (-1.*b*(-1.*r0 +b*cosbeta*costh/sinth)/sinth*sinbeta)/r0sq,                  //phi(0)
	 muin = (-1.*bsq*costh + 2.*r0sq*costh - 2.*b*r0*cosbeta*sinth)/2./r0sq,              //costh(0)=mu(0)
	 tprin=r0,                                                                            //t'(0)
	 rprin = bsq/(2.*r0) - r0,                                                            //r'(0)
	 muprin = (-1.*bsq*costh - 1.*b*r0*cosbeta*sinth)/r0sq,                               //mu'(0)
	 phprin=(b*(r0 - 2.*b*cosbeta*costh/sinth)/sinth*sinbeta)/r0sq;                       //phi'(0)
//initial perpendicular vector to a geodesic
doub tain=0.,                                                                             //t(0)
	 rain=-(b*cosbeta)/r0,                                                                //r(0)
	 muain = (2.*b*r0*cosbeta*costh - (bsq-2.*r0sq+2.*bsq*cosbetasq)*sinth)/(2.*r0sq*r0), //mu(0)
	 pain=(b*costh/sinth*(r0 - 2.*b*cosbeta*costh/sinth)/sinth*sinbeta)/r0sq*r0;          //phi(0)

t=0.;          //set initial affine parameters at zero
doub t1 = 2.1, //maximum affine parameter t1=2.1, which is never reached (even by geodesics going many times around the BH)
	 tmax;     //actual maximum affine parameter on a geodesic
//setting up primitive ODE variables
y[0] = tin;    
y[1]= phiin; 
y[2]=rprin; 
y[3]=muprin; 
y[4]=rin; 
y[5]=muin; 
y[6]=tprin;
y[7]=phprin;
y[8]=tain;
y[9]=rain;
y[10]=muain;
y[11]=pain;

solvegeodesic(t, y, der, NULL); //invoke geodesic equations once without a solver
tpast=t;
lam[0]=0.;                      //initial affine parameter
coox[0][stN]=y[0];              //coox array elements = primitive variables
coox[3][stN]=y[1];
coox[1][stN]=y[4];
coox[2][stN]=y[5];
coox[4][stN]=y[6];
coox[7][stN]=y[7];
coox[5][stN]=y[2];
coox[6][stN]=-y[3]/sqrt(1-y[5]*y[5]);
for(p=8;p<12;p++)
	coox[p][stN]=y[p];
for(p=0;p<12;p++){             //set global ppy variable = coox variable
	ppy[currth].cooxx[p][stN]=coox[p][stN];
}
ppy[currth].lamx[stN]=lam[stN];

//solver has
//y : t, phi, rp ,mup, r, mu ,tpr,phpr
//f : t',phi',r'',mu'',r',mu',t'',ph''
doub h=step;                   //initialize step of geodesic solver
//main integration cycle
while (t < t1){                //main geodesic integration loop
	stN++;                     //1 geodesic point per main loop
	int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t1, &h, y);//evolve geodesic solution with a chosen step
	lam[stN]=t;                //affine parameter at a next point
	coox[0][stN]=y[0];         //geodesic properties
	coox[3][stN]=y[1];
	coox[1][stN]=y[4];
	coox[2][stN]=y[5];
	coox[4][stN]=y[6];
	coox[7][stN]=y[7];
	coox[5][stN]=y[2];
	coox[6][stN]=-y[3]/sqrt(1-y[5]*y[5]);
	for(p=8;p<12;p++)
		coox[p][stN]=y[p];
	for(p=0;p<12;p++){         //set global ppy variable = coox variable
		ppy[currth].cooxx[p][stN]=coox[p][stN];
	}
	ppy[currth].lamx[stN]=lam[stN];
	if(tpast==t){              //check if step size is effectively zero - only happens if computing with reduced precision (float as opposed to double)
		printf("Geodesic cannot be computed due to low precision! Switch to double...\n Exiting");
		exit(-1);
	};
	tpast=t; 
	if(h>step)                 //if next suggested by evolve_apply step size is too large, then reduce it
		h=step;
	if ((status != GSL_SUCCESS) || (y[4]<(1.+ss)*rg)|| (y[4]>r0) || (stN>maxco))  //check if geodesic falls onto BH, exits the domain, or if too many points on geodesic are already computed 
		break;
}
if(stN>maxco){
	printf("Too many points requested at t=%e\n Exiting",t);
	exit(-1);
}
tmax=lam[stN];                //maximum affine parameter on a geodesic
txx=r0/3.;                    //auxiliary affine parameter at t=0.
ppy[currth].llmin=-log(txx+sqrt(1+txx*txx));//sinh^(-1) transformation at t=0. - correspondent minimum log-parameter
txx=r0*(tmax-1.)/3.;          //auxiliary affine parameter at t=tmax
ppy[currth].llmax=log(txx+sqrt(1+txx*txx));//sinh^(-1) transformation - correspondent maximum log-parameter
ppy[currth].indx=stN;         //number of points on geodesic written into global variable

gsl_odeiv_evolve_free (e);    //free memory
gsl_odeiv_control_free (c);
gsl_odeiv_step_free (s);