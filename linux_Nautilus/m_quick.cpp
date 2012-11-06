{doub xaccur,xaccurr,xfact,xss,xsnxy,xstep,xsstep,xIint,xIang,xans;int ind,fmin,fmax;string qadd="";
xaccur=3e-3; //+1. accuracy of geodesics computation - linear interpolation
xaccurr=1e-2;//+2. accuracy of radiative transfer integration
xfact=1.0;   //+3. size of the integration region
xss=1e-2;    //+4. fractional distance from BH horizon to start integration
xsnxy=111;   //+5. grid N in the picture plane
xstep=0.1;   //+6. step control of geodesic computation
xsstep=-0.09;//+7. step control of radiative transfer computation
xIint=2e-9;  //+8. initial intensity along each ray
xIang=1.;    //+9. initial polarized phases along each ray

//xsnxy=21;//for quick assessement
accur=xaccur;accurr=xaccurr;fact=xfact;ss=xss;int snxy=xsnxy;doub step=xstep,sstep=xsstep,Iint=xIint,Iang=xIang,devsq=0.;
cas=atoi(descr);//Bpo-=0.1*cas;
	start=std::clock();iswrite=true;
kmin=3;kmax=13;fdiff=60;fmin=6850;fmax=9850;sep=50;
kmin=2;kmax=13;
//for tests with intensity
for(kk=kmin;kk<=kmax;kk++){xtotin[kk]=0.;xLPo[kk]=0.;xCP[kk]=0.;xEVPA[kk]=0.;avrate=0.;avTpTe=0.;avTe6=0.;};

	if(sp==0){rhonor=1030713.25; heat=0.42107;  th=0.7328;fdiff=60;}//4.2802
    if(sp==1){rhonor= 988988.82; heat=0.37012;  th=1.8407;fdiff=60;}//2.9363
	if(sp==2){rhonor= 817113.52; heat=0.37239;  th=2.0159;fdiff=60;}//3.0837
	if(sp==3){rhonor= 656228.94; heat=0.39849;  th=2.2074;fdiff=60;}//4.1912
	if(sp==4){rhonor= 397281.19; heat=0.41343;  th=2.1439;fdiff=60;}//3.0059

	//if(sp==1){rhonor= 767789.57; heat=0.42077;  th=1.9056;fdiff=60;}//2.4738//update Mar2012 - doesn't pass convergence tests...
	
	if(sp==0){rhonor=1042075.25; heat=0.41929;  th=0.7231;fdiff=0;}//4.3485
	if(sp==1){rhonor= 872253.56; heat=0.39804;  th=1.8726;dense=1.;rhonor*=dense;fdiff=0;}//1.9967
	if(sp==2){rhonor= 884119.22; heat=0.35708;  th=2.1054;fdiff=0;}//2.4761
	if(sp==3){rhonor= 651827.84; heat=0.40215;  th=2.2109;fdiff=0;}//3.5874
	if(sp==4){rhonor= 385700.77; heat=0.41720;  th=2.1351;fdiff=0;}//3.5910 
	
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
case 17: fmin=4559;fmax=5559;sp=2;rhonor=1.e+7; heat=0.36;th=1.8407;fdiff=0;break;//ADAF simulations
case 18: fmin=6434;fmax=7034;sp=0;rhonor=1.05e+5; heat=0.17;th=1.9407;fdiff=0;break;//Jon simulations
case 19: fmin=6434;fmax=7034;sp=0;rhonor=1.05e+5; heat=0.17;th=2.9;fdiff=0;break;//Jon simulations
case 20: fmin=6434;fmax=7034;sp=0;rhonor=1.05e+5; heat=0.17;th=1.59;fdiff=0;break;//Jon simulations
case 21: fmin=6434;fmax=7034;sp=0;rhonor=5.e+4; heat=0.37;th=2.1;fdiff=0;break;//Jon simulations
case 22: fmin=5534;fmax=7034;sp=0;rhonor=144694.13267; heat=0.17339;th=2.4059;fdiff=0;thlimit=0.1;break;//Jon simulations with B^2/\rho cut-off
case 23: fmin=5534;fmax=7034;sp=0;rhonor=142642.42683; heat=0.16737;th=2.6367;thlimit=0.1;fdiff=0;isBcut=false;break;//Jon simulations, no B^2/\rho cut-off, thlimit=0.1
}printf("Bpo=%.3f, fdiff=%d\n",Bpo,fdiff);

ind=co;sep=(fmax-fmin)/(ind-1);
	tth=th;xth=th;
	for(fnum=fmin;fnum<=fmax;fnum+=sep){init(sp,fmin,fmax,sep);
	#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) shared(ittot)
	#include "intensity.cpp"
	for(kk=kmin;kk<=kmax;kk++){xtotin[kk]+=totin[kk]/ind;
		xLPo[kk]+=LPo[kk]/ind;
		xCP[kk]+=CP[kk]/ind;
		xEVPA[kk]+=ang[kk]/ind;
	};
	avrate+=rate/ind;avTpTe+=TpTe/ind;avTe6+=Te6/ind;
	};
	ans=( std::clock() - start) / (doub)CLOCKS_PER_SEC;printf ("Time = %.2f s; th=%.3f; heat=%.3f\n", ans,th,heat);


stringstream ytr;ytr<<(int)100*a<<"in"<<ind;string stra = ytr.str();
FILE * pFile; pFile = fopen ((dir+"quicka"+stra+qadd+".dat").c_str(),"a");

for(kk=kmin;kk<=kmax;kk++){
printf("avg at f=%.1f; I=%.3fJy LP=%.2f%% CP=%.3f%% EVPA=%.2f \n",sftab[kk][0], xtotin[kk],xLPo[kk], xCP[kk],xEVPA[kk]);
fprintf(pFile,"%.1f %.3f %.2f %.3f %.2f %.3f %.3e %.2f %.3e\n",sftab[kk][0], xtotin[kk],xLPo[kk], xCP[kk], xEVPA[kk],1.*cas,avrate*year/Msun,avTpTe,avTe6);
};

fclose(pFile);

;}
//}