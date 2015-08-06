time_t t_b4_geodesics = clock();

int p,                     //for loop variable
	stN=0;                 //counter of points along geodesic
doub tpast,                //time at a previous step
	 y[12],                //main array for a geodesic solver
	 der[12],              //auxiliary array of derivates, invoked by ODE solver //RG:derivates=derivatives?
	 txx,                  //auxiliary affine parameter for sinh conversion (see below)
	 lam[maxco],           //affine parameter along geodesic
	 coox[12][maxco];      //coordinates, tangential vector, and perpendicular parallel-transported vector in each point along geodesic
gsl_odeiv_step *s;         //auxiliary ODE variables from GSL (GNU Scientific Library)
gsl_odeiv_control *c; 
gsl_odeiv_evolve *e;

int currth = omp_get_thread_num(); //OpenMP thread ID => geodesic number //RG BUT: multiple threads per geodesic?

//printf(YELLOW"[geoint.cpp]: "RESET"accur=%g\n",accur); exit(1);

const gsl_odeiv_step_type * T;     //ODE solver type

//DEFAULT
T = gsl_odeiv_step_rk2;            //use Runge-Kutta 2-nd order
//RG:
//T = gsl_odeiv_step_rkf45; // include <gsl/...>

s = gsl_odeiv_step_alloc (T, 12);  //initialize with 12 variables
c = gsl_odeiv_control_standard_new(0.0, accur, 1.0, 0.0);//relative accuracy "accur"
e = gsl_odeiv_evolve_alloc (12);
gsl_odeiv_system geodesic_solver_gsl = {solvegeodesic, NULL, 12, NULL};//prepare to solve ODE of 12 variables given by "solvegeodesic" function

//a - spin, r0 - distance from picture plane to BH, th - cosine of polar angle, b - impact parameter, 
//beta - angle in a picture plane counterclockwise from the direction to the north pole  
//BH rotates clockwise as viewed from above the north pole

rg=1+sqrt(1.-asq);             //BH horizon //RG: rg is a terrible variable name rg->rh ?
doub bsq=b*b,                  //auxiliary geometric variables based on a, r0, th, b //RG: b=sqrt(xg*xg+yg*yg), //impact parameter
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
	 rin = 0.5*bsq/r0 + r0,                                                               //r(0) - not precisely r0 because of curvature, but assumed b<<r0
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

//RG: t~[0,2] geodesics show kinks for t~2
doub t1 = 2.1, //maximum affine parameter t1=2.1, which is never reached (even by geodesics going many times around the BH)
//doub t1 = 2.5, //maximum affine parameter t1=2.1, which is never reached (even by geodesics going many times around the BH)
	 tmax;     //actual maximum affine parameter on a geodesic

//RG: mathematica notebook: y : t, phi,  rp, mup,  r, mu, tpr, \[Phi]pr
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

for(p=8;p<12;p++) {

  // INITIALIZE ARRAY TO ZERO, still cooxx array has wrong values after geodesic is abandoned
  // for (int dummy_idx=0; dummy_idx<stN; dummy_idx++) {
  //   coox[p][dummy_idx]=0.;
  // }

  //RG: CHECK WHY DID WE SET coox[p][stN] just above to something else?
  coox[p][stN]=y[p];
}

for(p=0;p<12;p++){             //set global ppy variable = coox variable
	ppy[currth].cooxx[p][stN]=coox[p][stN];
    //geodesics[ix][stN][p] = coox[p][stN]; // STORE GEODESIC INFO IN TRULY GLOBAL ARRAY
}
ppy[currth].lamx[stN]=lam[stN];

//solver has
//y : t, phi, rp ,mup, r, mu ,tpr,phpr
//f : t',phi',r'',mu'',r',mu',t'',ph''
doub h=step;                   //initialize step of geodesic solver


/**********************************
 * main geodesic integration loop *
 **********************************/
while (t < t1){                //RG: t~[0,2] geodesics show kinks for t~2

  // INITIALIZE ARRAY TO ZERO, still cooxx array has wrong values after geodesic is abandoned
  // for (int dummy_idx=0; dummy_idx<stN; dummy_idx++) {
  //   coox[p][dummy_idx]=0.;
  // }


	stN++;                     //1 geodesic point per main loop
	int status = gsl_odeiv_evolve_apply (e, c, s, &geodesic_solver_gsl, &t, t1, &h, y);//evolve geodesic solution with a chosen step


    // doub d_to_image_fudge = 1.01;
	// if ((status != GSL_SUCCESS) || (y[4]<(1.+ss)*rg)|| (y[4]>d_to_image_fudge*r0) || (stN>maxco)) { //check if geodesic falls onto BH, exits the domain, or if too many points on geodesic are already computed 
    //   // printf(YELLOW"[geoint.cpp]: "RED"status=%d GSL_SUCCESS=%d y[4]=%g ss=%g rg=%g a=%f\n"RESET,status,GSL_SUCCESS,y[4],ss,rg,a);
    //   // if (y[4]<(1.+ss)*rg) {
    //   //   printf(YELLOW"[geoint.cpp]: "RED"EXITING BECAUSE GEODESIC GOT CLOSE TO THE HORIZON\n"RESET);
    //   //     }
    //   if (stN==1) printf(YELLOW"[geoint.cpp]: "RED"status=%d GSL_SUCCESS=%d y[4]=%g ss=%g rg=%g a=%f\n"RESET,status,GSL_SUCCESS,y[4],ss,rg,a);
    //   printf(YELLOW"[geoint.cpp]: "RED"EXIT WITH stN=%d\n"RESET,stN);
    //   break;
    // }

    

    //RG: DEBUG geodesic seems to be closer to equator than output diagnostic suggests
    int BAD_GEODESIC=4950; // 19980;
    doub z_geod=y[4]*y[5]; // z = r * costh

    // if ( /* fabs(costh)>1e-3 && */ ix==BAD_GEODESIC) printf(YELLOW"[geoint.cpp]: "RED"UGH OH! GEODESIC MOVING OFF EQUATORIAL PLANE! z=%g costh=%g geodesic_id=ix=%d h=%g stdN=%d\n"RESET,z_geod,y[5],ix,h,stN); //stN=1907 last point
    //if (ix==BAD_GEODESIC) printf(YELLOW"[geoint.cpp]: "RED"UGH OH! GEODESIC MOVING OFF EQUATORIAL PLANE! z=%g geodesic_id=ix=%d\n"RESET,z_geod,ix); 



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

    //printf(YELLOW"[geoint.cpp]: "RESET"affine parameter t=%f r=y[4]=%f stN=%d",t,y[4],stN);

	if(tpast==t){              //check if step size is effectively zero - only happens if computing with reduced precision (float as opposed to double)
		printf(YELLOW"[geoint.cpp]: "RED"Geodesic cannot be computed due to low precision! Switch to double...\nEXITING...\n"RESET);
		exit(-1);
	};
	tpast=t; 

	if(h>step) {                 //if next suggested by evolve_apply step size is too large, then reduce it //RG: WHY IS IT TOO LARGE? WHAT IS MEANT HERE? WHY ARE WE DOING THIS?
    //if (step<0) printf(YELLOW"[geoint.cpp]: "RED"Negative step size in geoint: step=%g h=%g\n"RESET,step,h);
    // if (false)

      // if (ix==0) printf(YELLOW"[geoint.cpp]: "RED"TIME STEP TOO LARGE! affine parameter t=%f h=%f step=%f"RESET"\n",t,h,step);
        h=step;
      }

    doub d_to_image_fudge = 1.01;

	if ((status != GSL_SUCCESS) || (y[4]<(1.+ss)*rg)|| (y[4]>d_to_image_fudge*r0) || (stN>maxco)) { //check if geodesic falls onto BH, exits the domain, or if too many points on geodesic are already computed 
      // printf(YELLOW"[geoint.cpp]: "RED"status=%d GSL_SUCCESS=%d y[4]=%g ss=%g rg=%g a=%f\n"RESET,status,GSL_SUCCESS,y[4],ss,rg,a);
      // if (y[4]<(1.+ss)*rg) {
      //   printf(YELLOW"[geoint.cpp]: "RED"EXITING BECAUSE GEODESIC GOT CLOSE TO THE HORIZON\n"RESET);
      //     }
      break;
    }

} // while (t < t1){


if(stN>maxco){
	printf(YELLOW"[geoint.cpp]: "RED"Too many points requested at t=%e\nIncrease maxco in [global_variables.cpp] \nEXITING\n"RESET,t);
	exit(-1);
}
tmax=lam[stN];                //maximum affine parameter on a geodesic

// RG: give equation ref 
// 1 + 3*Sinh[llog]/r0 == t // solve for llog
// l <-> t
// see [checks_GRMHD_code.nb] section: "Correct Dsolve, definition of lmin, lmax, interpolation"
txx=r0/3.;                    //auxiliary affine parameter at t=0. // r0 distance to picture plane=maximum distance to BH of point on geodesic
ppy[currth].llmin=-log(txx+sqrt(1+txx*txx));//sinh^(-1) transformation at t=0. - correspondent minimum log-parameter
txx=r0*(tmax-1.)/3.;          //auxiliary affine parameter at t=tmax
ppy[currth].llmax=log(txx+sqrt(1+txx*txx));//sinh^(-1) transformation - correspondent maximum log-parameter
ppy[currth].indx=stN;         //number of points on geodesic written into global variable


/***************************/
/* RG:OUTPUT GEODESIC INFO */

// GEODESIC DIAGNOSTIC
stringstream geodesic_sstr;
FILE * pFile; 

geodesic_sstr<<"geodesics"<<(int)100*a<<"th"<<(int)floor(100*th+1e-6)<<"fn"<<fnum<<"geodesic";
string append_label;
append_label=geodesic_sstr.str();

// if (iix%geodesic_output_every_x==0 && iiy%geodesic_output_every_y==0) { // sample full 2d picture plane
// if (iix%geodesic_output_every_x==0 && iiy==nxy/2) { // sample y=const slice of picture plane
if (iiy%geodesic_output_every_x==0 && iix==nxy/2) { // sample x=const slice of picture plane

  int geodesic_label=ix;
  stringstream geodesic_sstr_append;
  geodesic_sstr_append<<append_label<<geodesic_label;
  string stra = geodesic_sstr_append.str();

  pFile=fopen ((dir+stra+".dat").c_str(),"a");

  //for (int geo_idx=0; geo_idx<=maxco; geo_idx++) { // maxco
  for (int geo_idx=0; geo_idx<=stN; geo_idx++) { // maxco
    for (int p=0; p<=7; p++) { // r:p-> costh: ph:
      fprintf(pFile,"%f ",ppy[currth].cooxx[p][geo_idx]);
    }
    // lamx[] vs lam[]?
    fprintf(pFile,"%f %d %d\n",ppy[currth].lamx[geo_idx],geodesic_label,geo_idx);
  }

  fclose(pFile);

} //for (int geodesic_label=0; geo_idx<=maxco; geo_idx+=1000) { // if (currth)


/********* END OF GEODESIC DIAGNOSTIC ************/


/***************************/


gsl_odeiv_evolve_free (e);    //free memory
gsl_odeiv_control_free (c);
gsl_odeiv_step_free (s);

t_geodesics += (clock() - t_b4_geodesics) / (doub)CLOCKS_PER_SEC;
