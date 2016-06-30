{//image of the emitting region
//calculations are performed with better accuracy compared to other m_xxx.cpp calculations

//RG: DEFAULT:
doub xaccur=3e-4,  //1. absolute accuracy of geodesics computation
 	 xaccurr=5e-4, //2. absolute accuracy of radiative transfer integration
//RG: 1e-6 too many points requested...

	 xfact=1.0,    //3. relative size of the integration region
     xss=3e-3,     //4. fractional distance from BH horizon to the sphere, where geodesic integration stops
     // xss=5e-1,     //4. fractional distance from BH horizon to the sphere, where geodesic integration stops
     xsnxy=nxy,  //5. number of points N along one side in the picture plane for N x N intensity grid
     xstep=1e-2,   //6. step size in geodesic computation
     // xstep=1e-3,   //6. step size in geodesic computation
     xsstep=-0.06, //7. step size in radiative transfer computation
     // xsstep=-0.01, //7. step size in radiative transfer computation
	 xIint=1e-10,  //8. initial intensity along each ray for radiative transfer
	 xIang=0.1;    //9. initial polarized phases along each ray for radiative transfer
doub step=xstep,   //local variables, which control radiative transfer
	 sstep=xsstep,
	 Iint=xIint,
	 Iang=xIang;
int ind,           //number of time frames for spectrum evaluation (to compute mean spectrum)
	fmin,          //minimum ID of fluid simulation snapshot
	fmax;          //maximum ID of fluid simulation snapshot

accur=xaccur;      //assigning values to global variables, which control radiative transfer
accurr=xaccurr;
fact=xfact;
ss=xss;

kmin=0;kmax=10;    //define a set of frequencies from "sftab" table
kstep=1;           //step in an array of frequencies

cas=atoi(descr);   //batch job ID = selection of model spin, snapshot IDs, density, heating constant, viewing angle + auxiliary parameters
start=clock(); 
iswrite=true;      //write results to disk

switch (cas){      //selection of a model (only few examples are shown)
	case 0: fmin=6850;fmax=9850; sp=1; rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;break;                                   //sample model without any changes to temperature in hot/tenuous regions
	case 7: fmin=6850;fmax=9850; sp=1; rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;fif="_all";break;                       //model with all radiative transfer effects on
	case 8: fmin=6850;fmax=9850; sp=1; rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;fif="_jVc0";fljVc=0.;break;             //V-mode emissivity is set to zero
	case 9: fmin=6850;fmax=9850; sp=1; rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;fif="_rQc0";flrQc=0.;break;             //Faraday conversion is set to zero
	case 10: fmin=6850;fmax=9850;sp=1; rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;fif="_rVc0";flrVc=0.;break;             //Faraday rotation is set to zero
	case 11: fmin=6850;fmax=9850;sp=1; rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;fif="_Bp15";Bpo-=0.5;break;             //change magnetic field extension slope - first synchronize with command line arguments!
	case 12: fmin=10000;fmax=18720;sp=0;rhonor=2.5e+8; heat=0.55;th=1.8407;fdiff=0;break;                                         //fast light approximation
	case 14: sp=0;fmin=12424;fmax=22424;rhonor=14802376.9520; heat=0.59894;th=2.356;thlimit=0.1;isBcut=false;fdiff=0;break;       //introducing non-zero thlimit, setting isBcut and isBred
	case 35: sp=0;rhonor=261385.84479; heat=0.15572;th=1.745/*10 deg from edge-on*/;thlimit=0.1;fdiff=30;isBcut=false;isBred=true;break;//changing viewing angle, almost edge-on
	case 36: sp=0;rhonor=261385.84479; heat=0.15572;th=2.967/*10 deg from face-on*/;thlimit=0.1;fdiff=30;isBcut=false;isBred=true;break;//changing viewing angle, almost face-on
	case 37: sp=0;rhonor=147780.66084; heat=0.16992;th=2.4840;dphi=2.*PI/3.;thlimit=0.1;fdiff=30;isBcut=false;isBred=false;break;//changing azimuthal camera angle - 2*PI/3
	case 38: sp=0;rhonor=147780.66084; heat=0.16992;th=2.4840;dphi=4.*PI/3.;thlimit=0.1;fdiff=30;isBcut=false;isBred=false;break;//changing azimuthal camera angle - 4*PI/3
    //RG:
      //case 112: fmin=6100;fmax=6230;kmin=0;kmax=0;sp=0;rhonor=147780.66084; heat=0.16992;th=1.745/*2.4840*/;dphi=4.*PI/3.;thlimit=0.0;fdiff=0;isBcut=false;isBred=false;break;//changing azimuthal camera angle - 4*PI/3
      //case 112: fmin=6100;fmax=6230;kmin=0;kmax=0;sp=0;rhonor=147780.66084; heat=0.16992;th=1.745/*2.4840*/;dphi=4.*PI/3.;thlimit=0.0;fdiff=0;isBcut=false;isBred=true;break;//changing azimuthal camera angle - 4*PI/3
      //case 112: fmin=6100;fmax=6230;kmin=0;kmax=0;sp=0;rhonor=147780.66084; heat=0.16992;th=1.745/*2.4840*/;dphi=4.*PI/3.;thlimit=0.0;fdiff=0;isBcut=true;isBred=false;break;//changing azimuthal camera angle - 4*PI/3
      //case 112: fmin=6100;fmax=6230;kmin=0;kmax=0;sp=0;rhonor=147780.66084; heat=0.16992;th=1.745/*2.4840*/;dphi=4.*PI/3.;thlimit=0.1;fdiff=0;isBcut=true;isBred=false;break;//changing azimuthal camera angle - 4*PI/3
      //case 112: fmin=6100;fmax=6200;kmin=0;kmax=0;sp=0;rhonor=147780.66084; heat=0.16992; th=0.5; dphi=4.*PI/3.;thlimit=0.1;fdiff=0;isBcut=true;isBred=false;break;//changing azimuthal camera angle - 4*PI/3

      // RG: USE THESE FILES FOR VARYING PARAMETERS (MOVIES, ETC) AND KEEP [m_imag.cpp] TIDY
      //  #include "vary_theta_slices.cpp"
      //#include "vary_r_slices.cpp"
      //#include "vary_magn.cpp"
      //#include "vary_thlimit.cpp"
      #include "models.cpp"
};
printf(YELLOW"[m_imag.cpp]: "RESET"Bpo=%.3f fdiff=%d th=%f\n",Bpo,fdiff,th);


/* SETTING UP MAIN LOOP OVER SNAPSHOTS */
//ind=21;                 //set a number of fluid simulation snapshots (an image is computed for each) (average spectra over ind nr of snapshots? see appendix in Shcherbakov,Penna & McKinney 2012)
//sep=(fmax-fmin)/(ind-1);//compute difference of IDs of two consecutive snapshots

//RG: catch fmin=fmax
if (sep==0) sep=1;

// WHY NOT JUST USE 
int fnum_step=10;

for(fnum=fmin;fnum<=fmax;fnum+=sep){//compute images, parallelization with OpenMP
  
  //RG:
  printf(YELLOW"[m_imag.cpp]: "RESET"fnum=%d,sep=%d\n",fnum,sep);
  
  init(sp,fmin,fmax,sep); //RG: sep not used in init()... remove!

  // OPTION TO CYCLE OVER CAMERA VIEWING ANGLE, INCLINATION, ETC FOR GIVEN CASE
  const string SCAN_THROUGH=""; // "view";
  doub SCAN_THROUGH_MAX=0.;
  if      (SCAN_THROUGH=="view") {
    SCAN_THROUGH_MAX=2.*PI; // [0,2pi] in dphi
    dphi=0.;
  }
  else if (SCAN_THROUGH=="inclination") {
    SCAN_THROUGH_MAX=PI;    // [0,pi] in th
    th=0.;                  // MIGHT WANT TO AVOID POLES... :-s
  }
  else
    SCAN_THROUGH_MAX=0;    // DEACTIVATE


  for (doub scan_through_value=0;scan_through_value<=SCAN_THROUGH_MAX;scan_through_value+=max(SCAN_THROUGH_MAX/100.,1e-9)){

    if (SCAN_THROUGH=="view") dphi=0.+scan_through_value;
    
    if (SCAN_THROUGH!="") {
      printf(YELLOW"[m_imag.cpp]: "RESET"SCANNING THROUGH %s, currently at %f...\n",SCAN_THROUGH.c_str(),scan_through_value);
    }

    // RG: IS DYNAMIC SCHEDULING REALLY NEEDED?
#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) shared(ittot)
#include "imaging.cpp"	

  } // SCAN_THROUGH LOOP

}
ans=(clock() - start) / (doub)CLOCKS_PER_SEC;printf ("Time = %.2f s; th=%.3f; heat=%.3f\n", ans,th,heat);

}
