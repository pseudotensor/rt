{//quick computation of spectrum for given sp, heat, rhonor, th
// doub xaccur=3e-3,  //1. absolute accuracy of geodesics computation
// 	 xaccurr=1e-3, //2. absolute accuracy of radiative transfer integration
// xfact=1.0,    //3. relative size of the integration region
// xss=1e-2,     //4. fractional distance from BH horizon to the sphere, where geodesic integration stops
	 // xsnxy=151,    //5. number of points N along one side in the picture plane for N x N intensity grid
	 // xstep=0.1,    //6. step size in geodesic computation
	 // xsstep=-0.09, //7. step size in radiative transfer computation
	 // xIint=2e-9,   //8. initial intensity along each ray for radiative transfer
	 // xIang=1.;     //9. initial polarized phases along each ray for radiative transfer

// as in original ASTRORAY m_imag
doub xaccur=3e-4,  //1. absolute accuracy of geodesics computation
	 xaccurr=5e-4, //2. absolute accuracy of radiative transfer integration
     xfact=1.0,    //3. relative size of the integration region
	 xss=3e-3,     //4. fractional distance from BH horizon to the sphere, where geodesic integration stops

  // xsnxy=101,    //5. number of points N along one side in the picture plane for N x N intensity grid
	 xsnxy=nxy,    //5. number of points N along one side in the picture plane for N x N intensity grid
	 xstep=1e-2,    //6. step size in geodesic computation
	 xsstep=-0.06, //7. step size in radiative transfer computation
	 xIint=1e-10,   //8. initial intensity along each ray for radiative transfer
	 xIang=0.1;     //9. initial polarized phases along each ray for radiative transfer
doub step=xstep,   //local variables, which control radiative transfer
	 sstep=xsstep,
	 Iint=xIint,
	 Iang=xIang;
doub avrate,       //average accretion rate
	 avTpTe,       //average ratio of Tp/Te at 6M radius
	 avTe6;        //average electron temperature Te at 6M radius
//int snxy=xsnxy,    //global variable correspondent to xsnxy
int snxy=nxy,    //global variable correspondent to xsnxy
    ind,           //number of time frames for spectrum evaluation // RG: average over ind snapshots?
	fmin,          //minimum ID of fluid simulation snapshot
	fmax,          //maximum ID of fluid simulation snapshot
	sep;           //ID difference between consecutive considered fluid simulation snapshots

accur=xaccur;      //assigning values to global variables, which control radiative transfer
accurr=xaccurr;
fact=xfact;
ss=xss;

cas=atoi(descr);   //batch job ID = selection of model spin, snapshot IDs, density, heating constant, viewing angle + auxiliary parameters
start=clock();
iswrite=true;      //write results to disk
kmin=0;kmax=13;    //selection of minimum and maximum frequencies from "sftab" array
fmin=6850;fmax=9850;//default minimum and maximum fluid simulation IDs
sep=50;            //default ID difference between consecutively taken fluid simulation snapshots
fdiff=30;          //default d(ID) so that light rays propagate as simulation evolves from fnum-fdiff to fnum+fdiff IDs

                   //initializing with zeros
for(kk=kmin;kk<=kmax;kk++){
	xtotin[kk]=0.;
	xLPo[kk]=0.;
	xCP[kk]=0.;
	xEVPA[kk]=0.;
	avrate=0.;
	avTpTe=0.;
	avTe6=0.;
};

switch (cas){      //selection of a model (only few examples are shown)
	case 0: fmin=6850;fmax=9850; sp=1; rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;break;                                   //sample model without any changes to temperature in hot/tenuous regions
	case 7: fmin=6850;fmax=9850; sp=1; rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;fif="_all";break;                       //model with all radiative transfer effects on
	case 8: fmin=6850;fmax=9850; sp=1; rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;fif="_jVc0";fljVc=0.;break;             //V-mode emissivity is set to zero
	case 9: fmin=6850;fmax=9850; sp=1; rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;fif="_rQc0";flrQc=0.;break;             //Faraday conversion is set to zero
	case 10: fmin=6850;fmax=9850;sp=1; rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;fif="_rVc0";flrVc=0.;break;             //Faraday rotation is set to zero
	case 11: fmin=6850;fmax=9850;sp=1; rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;fif="_Bp15";Bpo-=0.5;break;             //change magnetic field extension slope - first synchronize with command line arguments!
	case 12: fmin=10000;fmax=18720;sp=0;rhonor=2.5e+8; heat=0.55;th=1.8407;fdiff=0;break; //fast light approximation // RG: fdifff=0 and fmin!=fmax does not make sense to me...
      // case 112: fmin=7033;fmax=7033;sp=0;rhonor=2.5e+8; heat=0.55;th=1.8407;fdiff=0;break;
    // case 112: fmin=6900;fmax=7016;sp=0; rhonor=14802376.9520; heat=0.59894;th=2.356;thlimit=0.1;isBcut=false;fdiff=0;kmin=7;kmax=7;break;
// RG: I introduce my own cases and use 3 digit cases // look below at include commands 

	case 14: sp=0;fmin=12424;fmax=22424;rhonor=14802376.9520; heat=0.59894;th=2.356;thlimit=0.1;isBcut=false;fdiff=0;break;       //introducing non-zero thlimit, setting isBcut and isBred
	case 35: sp=0;rhonor=261385.84479; heat=0.15572;th=1.745/*10 deg from edge-on*/;thlimit=0.1;fdiff=30;isBcut=false;isBred=true;break;//changing viewing angle, almost edge-on
	case 36: sp=0;rhonor=261385.84479; heat=0.15572;th=2.967/*10 deg from face-on*/;thlimit=0.1;fdiff=30;isBcut=false;isBred=true;break;//changing viewing angle, almost face-on
	case 37: sp=0;rhonor=147780.66084; heat=0.16992;th=2.4840;dphi=2.*PI/3.;thlimit=0.1;fdiff=30;isBcut=false;isBred=false;break;//changing azimuthal camera angle - 2*PI/3
	case 38: sp=0;rhonor=147780.66084; heat=0.16992;th=2.4840;dphi=4.*PI/3.;thlimit=0.1;fdiff=30;isBcut=false;isBred=false;break;//changing azimuthal camera angle - 4*PI/3

      #include "lightup_jet.cpp"
};


//sep=1;                                           //take each fluid simulation snapshot within the range
ind=(fmax-fmin)/sep+1;                           //number of fluid simulation snapshots

printf(YELLOW"[m_quick.cpp]: "RESET"cas=%d,th=%f,rhonor=%f,heat=%f,ind=%d,fmin=%d,fmax=%d,sep=%d,Bpo=%.3f, fdiff=%d\n",cas,th,rhonor,heat,ind,fmin,fmax,sep,Bpo,fdiff);



for(fnum=fmin;fnum<=fmax;fnum+=sep){             //cycle over fluid simulation snapshots - computing average spectrum

    cout << YELLOW"[m_quick.cpp]:"RESET" CALLING init()...\n";
	init(sp,fmin,fmax,sep);

    //cout << YELLOW"[m_quick.cpp]:"RESET"include intensity.cpp\n";
	#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) shared(ittot)
	#include "intensity.cpp"
	for(kk=kmin;kk<=kmax;kk++){
		xtotin[kk]+=totin[kk]/ind;
		xLPo[kk]+=LPo[kk]/ind;
		xCP[kk]+=CP[kk]/ind;
		xEVPA[kk]+=EVPA[kk]/ind;
	};
	avrate+=rate/ind;                            //computing averages of other quantities
	avTpTe+=TpTe/ind;
	avTe6+=Te6/ind;
};
ans=(clock() - start) / (doub)CLOCKS_PER_SEC;
// printf(YELLOW"[m_quick.cpp]:"RESET" Time = %.2f s; th=%.3f; heat=%.3f\n", ans,th,heat);

 stringstream ytr;                               //writing average spectrum into "quicka" file, specifying "ind"
 stringstream rt_pars;
ytr<<(int)100*a<<"in"<<ind<<"case"<<cas;
rt_pars<<"th"<<(float)th<<"rhonor"<<(float)rhonor<<"Cheat"<<heat;
string stra = ytr.str();
FILE * pFile; 

//RG: 
//cout << YELLOW"[m_quick.cpp]:"RESET" FILE:"+dir+"quicka"+stra+fif+".dat" << "\n";

 pFile = fopen ((dir+"quicka"+stra+rt_pars.str()+fif+".dat").c_str(),"a");


for(kk=kmin;kk<=kmax;kk++){                     //actual writing into "quicka" file
  	printf(YELLOW"[m_quick.cpp]:"RESET"  avg at f=%.1f; I=%.3fJy LP=%.2f%% CP=%.3f%% EVPA=%.2f \n",sftab[kk][0], xtotin[kk],xLPo[kk], xCP[kk],xEVPA[kk]);
	fprintf(pFile,"%.1f %.3f %.2f %.3f %.2f %.3f %.3e %.2f %.3e\n",sftab[kk][0], xtotin[kk],xLPo[kk], xCP[kk], xEVPA[kk],1.*cas,avrate*year/Msun,avTpTe,avTe6);
};
fclose(pFile);
};
