{//image of the emitting region
//calculations are performed with better accuracy compared to other m_xxx.cpp calculations
doub xaccur=3e-4,  //1. absolute accuracy of geodesics computation
	 xaccurr=5e-4, //2. absolute accuracy of radiative transfer integration
	 xfact=1.0,    //3. relative size of the integration region
	 xss=3e-3,     //4. fractional distance from BH horizon to the sphere, where geodesic integration stops
	 //xsnxy=101,  //5. number of points N along one side in the picture plane for N x N intensity grid
	 xstep=1e-2,   //6. step size in geodesic computation
	 xsstep=-0.06, //7. step size in radiative transfer computation
	 xIint=1e-10,  //8. initial intensity along each ray for radiative transfer
	 xIang=0.1;    //9. initial polarized phases along each ray for radiative transfer
doub step=xstep,   //local variables, which control radiative transfer
	 sstep=xsstep,
	 Iint=xIint,
	 Iang=xIang;
int ind,           //number of time frames for spectrum evaluation
	fmin,          //minimum ID of fluid simulation snapshot
	fmax;          //maximum ID of fluid simulation snapshot
string qadd="";    //output filename modifier, helps to distinguish cases

accur=xaccur;      //assigning values to global variables, which control radiative transfer
accurr=xaccurr;
fact=xfact;
ss=xss;

kmin=2;kmax=10;    //define a set of frequencies from "sftab" table
kstep=2;           //step in an array of frequencies

cas=atoi(descr);   //batch job ID = selection of model spin, snapshot IDs, density, heating constant, viewing angle + auxiliary parameters
start=clock(); 
iswrite=true;      //write results to disk

switch (cas){      //selection of a model (only few examples are shown)
	case 0: fmin=6850;fmax=9850; sp=1; rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;break;                                   //sample model without any changes to temperature in hot/tenuous regions
	case 7: fmin=6850;fmax=9850; sp=1; rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;qadd="_all";break;                       //model with all radiative transfer effects on
	case 8: fmin=6850;fmax=9850; sp=1; rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;qadd="_jVc0";fljVc=0.;break;             //V-mode emissivity is set to zero
	case 9: fmin=6850;fmax=9850; sp=1; rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;qadd="_rQc0";flrQc=0.;break;             //Faraday conversion is set to zero
	case 10: fmin=6850;fmax=9850;sp=1; rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;qadd="_rVc0";flrVc=0.;break;             //Faraday rotation is set to zero
	case 11: fmin=6850;fmax=9850;sp=1; rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;qadd="_Bp15";Bpo-=0.5;break;             //change magnetic field extension slope - first synchronize with command line arguments!
	case 12: fmin=10000;fmax=18720;sp=0;rhonor=2.5e+8; heat=0.55;th=1.8407;fdiff=0;break;                                         //fast light approximation
	case 14: sp=0;fmin=12424;fmax=22424;rhonor=14802376.9520; heat=0.59894;th=2.356;thlimit=0.1;isBcut=false;fdiff=0;break;       //introducing non-zero thlimit, setting isBcut and isBred
	case 35: sp=0;rhonor=261385.84479; heat=0.15572;th=1.745/*10 deg from edge-on*/;thlimit=0.1;fdiff=30;isBcut=false;isBred=true;break;//changing viewing angle, almost edge-on
	case 36: sp=0;rhonor=261385.84479; heat=0.15572;th=2.967/*10 deg from face-on*/;thlimit=0.1;fdiff=30;isBcut=false;isBred=true;break;//changing viewing angle, almost face-on
	case 37: sp=0;rhonor=147780.66084; heat=0.16992;th=2.4840;dphi=2.*PI/3.;thlimit=0.1;fdiff=30;isBcut=false;isBred=false;break;//changing azimuthal camera angle - 2*PI/3
	case 38: sp=0;rhonor=147780.66084; heat=0.16992;th=2.4840;dphi=4.*PI/3.;thlimit=0.1;fdiff=30;isBcut=false;isBred=false;break;//changing azimuthal camera angle - 4*PI/3
};
printf("Bpo=%.3f, fdiff=%d\n th=%f\n",Bpo,fdiff,th);

ind=21;                 //set a number of fluid simulation snapshots (an image is computed for each)
sep=(fmax-fmin)/(ind-1);//compute difference of IDs of two consecutive snapshots

for(fnum=fmin;fnum<=fmax;fnum+=sep){//compute images, parallelization with OpenMP
	init(sp,fmin,fmax,sep);
	#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) shared(ittot)
	#include "imaging.cpp"	
}
ans=(clock() - start) / (doub)CLOCKS_PER_SEC;printf ("Time = %.2f s; th=%.3f; heat=%.3f\n", ans,th,heat);
}