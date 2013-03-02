{doub xaccur,xaccurr,xfact,xss,xsnxy,xstep,xsstep,xIint,xIang,xans;int m,ind,fmin,fmax;string qadd="";
//better accuracy than for the rest of computations
xaccur=3e-4; //+1. accuracy of geodesics computation
xaccurr=5e-4;//+2. accuracy of radiative transfer integration
xfact=1.0;   //3. size of the integration region
xss=3e-3;    //+4. fractional distance from BH horizon to start integration
//xsnxy=101;    //+5. grid N in the picture plane//
xstep=1e-2;   //+6. step control of geodesic computation
xsstep=-0.06;//+7. step control of radiative transfer computation
xIint=1e-10;  //+8. initial intensity along each ray
xIang=0.1;    //+9. initial polarized phases along each ray
accur=xaccur;accurr=xaccurr;fact=xfact;ss=xss;doub step=xstep,sstep=xsstep,Iint=xIint,Iang=xIang,devsq=0.;
//int nxy=xsnxy;//should be no more than 401
//! using h
	start=std::clock();iswrite=true;
kmin=4;kmax=10;//for imaging 4-> 87GHz, 7 ->230GHz; 10 -> 857GHz
kmin=4;kmax=10;kstep=3;
// 2 0 2.0 4  -- are command lines

//fmin=6850;fmax=9850;sep=150;//ind=0..20
cas=atoi(descr);

start=std::clock();
	if(sp==0){rhonor=1030713.25; heat=0.42107;  th=0.7328;fdiff=60;}//4.2802
    if(sp==1){rhonor= 988988.82; heat=0.37012;  th=1.8407;fdiff=60;}//2.9363
	if(sp==2){rhonor= 817113.52; heat=0.37239;  th=2.0159;fdiff=60;}//3.0837
	if(sp==3){rhonor= 656228.94; heat=0.39849;  th=2.2074;fdiff=60;}//4.1912
	if(sp==4){rhonor= 397281.19; heat=0.41343;  th=2.1439;fdiff=60;}//3.0059

switch (cas){
case 0: fmin=6850;fmax=9850;if(sp==1){rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;};break;//2.4738
case 1: fmin=6850;fmax=7350;if(sp==1){rhonor=756365.55; heat=0.38934;th=1.9202;fdiff=60;};break;//4.5442
case 2: fmin=7350;fmax=7850;if(sp==1){rhonor=743209.25; heat=0.40507;th=1.8716;fdiff=60;};break;//5.1106
case 3: fmin=7850;fmax=8350;if(sp==1){rhonor=838898.40; heat=0.37302;th=1.8609;fdiff=60;};break;//3.3787
case 4: fmin=8350;fmax=8850;if(sp==1){rhonor=970599.66; heat=0.36147;th=1.8431;fdiff=60;};break;//3.1890
case 5: fmin=8850;fmax=9350;if(sp==1){rhonor=1148284.41;heat=0.37420;th=1.8871;fdiff=60;};break;//2.3095
case 6: fmin=9350;fmax=9850;if(sp==1){rhonor=1309041.71;heat=0.38853;th=1.8080;fdiff=60;};break;//3.0481
	//switching on/off various effects
case 7: fmin=6850;fmax=9850;if(sp==1){rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;};qadd="_all";break;
case 8: fmin=6850;fmax=9850;if(sp==1){rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;};qadd="_jVc0";fljVc=0.;break;
case 9: fmin=6850;fmax=9850;if(sp==1){rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;};qadd="_rQc0";flrQc=0.;break;
case 10: fmin=6850;fmax=9850;if(sp==1){rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;};qadd="_rVc0";flrVc=0.;break;
case 11: fmin=6850;fmax=9850;if(sp==1){rhonor=988988.82; heat=0.37012;th=1.8407;fdiff=60;};qadd="_Bp15";Bpo-=0.5;break;
case 12: fmin=10000;fmax=18720;sp=0;rhonor=2.5e+8; heat=0.55;th=1.8407;fdiff=0;break;//new simulations
case 13: fmin=8999;fmax=10299;sp=2;rhonor=6.e+7; heat=0.55;th=1.7;fdiff=0;break;//new simulations
case 14: fmin=4000;fmax=5500;sp=3;rhonor=2.4e+7; heat=0.55;th=1.7;fdiff=0;break;//new simulations
case 15: fmin=2700;fmax=2750;sp=4;rhonor=2.4e+7; heat=0.55;th=1.7;fdiff=0;break;//new simulations
case 16: fmin=10095;fmax=11747;sp=0;rhonor=1.e+7; heat=0.32;th=1.8407;fdiff=0;break;//ADAF simulations
case 17: fmin=4459;fmax=5559;sp=2;rhonor=1.e+7; heat=0.36;th=1.8407;fdiff=0;break;//ADAF simulations
case 18: fmin=6434;fmax=7034;sp=0;rhonor=1.05e+5; heat=0.17;th=1.9407;fdiff=0;break;//Jon simulations
case 19: fmin=6434;fmax=7034;sp=0;rhonor=1.05e+5; heat=0.17;th=2.9;fdiff=0;break;//Jon simulations
case 20: fmin=6434;fmax=7034;sp=0;rhonor=1.05e+5; heat=0.17;th=1.59;fdiff=0;break;//Jon simulations
case 21: fmin=6434;fmax=7034;sp=0;rhonor=5.e+4; heat=0.37;th=2.1;fdiff=0;break;//Jon simulations
case 22: fmin=5534;fmax=7034;sp=0;rhonor=144694.13267; heat=0.17339;th=2.4059;fdiff=0;break;//Jon simulations with B^2/\rho cut-off
case 23: fmin=5534;fmax=7034;sp=0;rhonor=142642.42683; heat=0.16737;th=2.6367;thlimit=0.1;fdiff=0;isBcut=false;break;//Jon simulations, no B^2/\rho cut-off, thlimit=0.1
case 24: fmin=5534;fmax=7034;sp=0;rhonor=274614.53874; heat=0.14422;th=2.4483;thlimit=0.1;fdiff=0;isBcut=true;break;//Jon simulations, both B^2/\rho cut-off and thlimit=0.1
case 28: fmin=5534;fmax=7034;sp=0;rhonor=288331.85590; heat=0.14320;th=2.4529;thlimit=0.1;fdiff=0;isBcut=false;isBred=true;break;//Jon simulations, thlimit=0.1 and rho reduction; no B^2/\rho cut-off
case 29: sp=0;rhonor=261385.84479; heat=0.15572;th=2.3768;thlimit=0.1;fdiff=30;isBcut=false;isBred=true;break;//Jon simulations, thlimit=0.1 and rho reduction; no B^2/\rho cut-off
case 30: sp=0;rhonor=147780.66084; heat=0.16992;th=2.4840;thlimit=0.1;fdiff=30;isBcut=false;isBred=false;break;//Jon simulations, thlimit=0.1; no rho reduction and no B^2/\rho cut-off
case 31: sp=0;rhonor=147780.66084; heat=0.16992;th=1.745/*10 deg from edge-on*/;thlimit=0.1;fdiff=30;isBcut=false;isBred=false;break;//Jon simulations, thlimit=0.1; no rho reduction and no B^2/\rho cut-off
case 32: sp=0;rhonor=147780.66084; heat=0.16992;th=2.967/*10 deg from face-on*/;thlimit=0.1;fdiff=30;isBcut=false;isBred=false;break;//Jon simulations, thlimit=0.1; no rho reduction and no B^2/\rho cut-off
case 33: sp=0;rhonor=261385.84479; heat=0.15572;th=2.3768;thlimit=0.1;fdiff=0;isBcut=false;isBred=true;break;//Jon simulations, thlimit=0.1 and rho reduction; no B^2/\rho cut-off
case 34: sp=0;rhonor=144897.83871; heat=0.17172;th=2.3897;thlimit=0.1;fdiff=0;isBcut=false;isBred=false;break;//Jon simulations, thlimit=0.1; no rho reduction and no B^2/\rho cut-off
case 35: sp=0;rhonor=261385.84479; heat=0.15572;th=1.745/*10 deg from edge-on*/;thlimit=0.1;fdiff=30;isBcut=false;isBred=true;break;//Jon simulations, thlimit=0.1 and rho reduction; no B^2/\rho cut-off
case 36: sp=0;rhonor=261385.84479; heat=0.15572;th=2.967/*10 deg from face-on*/;thlimit=0.1;fdiff=30;isBcut=false;isBred=true;break;//Jon simulations, thlimit=0.1 and rho reduction; no B^2/\rho cut-off
};printf("Bpo=%.3f, fdiff=%d at frequency f=%.3fGHz\n",Bpo,fdiff,sftab[kmax][0]);
printf("th=%f\n",th);tth=th;xth=th;

fmin=co;fmax=mco;sep=1;
//ind=co;sep=(fmax-fmin)/(ind-1);

//for(fnum=0;m<ind;m++){//comment out for ind given by PBS_ARRAYID
for(fnum=fmin;fnum<=fmax;fnum+=sep){
//fnum=fmin+m*sep;
	init(sp,fmin,fmax,sep);
	#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) shared(ittot)
	#include "imaging.cpp"	
}//comment out for ind given by PBS_ARRAYID
	ans=( std::clock() - start) / (doub)CLOCKS_PER_SEC;printf ("Time = %.2f s; th=%.3f; heat=%.3f\n", ans,th,heat);
}