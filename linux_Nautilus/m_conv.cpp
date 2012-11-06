{doub xaccur,xaccurr,xfact,xss,xsnxy,xstep,xsstep,xIint,xIang,xans;int testN, ind;
xaccur=3e-3; //+1. accuracy of geodesics computation
xaccurr=1e-2;//+2. accuracy of radiative transfer integration
xfact=1.0;   //3. size of the integration region
xss=1e-2;    //+4. fractional distance from BH horizon to start integration
xsnxy=111;    //+5. grid N in the picture plane
xstep=0.1;   //+6. step control of geodesic computation
xsstep=-0.09;//+7. step control of radiative transfer computation
xIint=2e-9;  //+8. initial intensity along each ray
xIang=1.;    //+9. initial polarized phases along each ray

//xsnxy=21;//for quick tests
//1 21 2.0 4 - command line

accur=xaccur;accurr=xaccurr;fact=xfact;ss=xss;int snxy=xsnxy;doub step=xstep,sstep=xsstep,Iint=xIint,Iang=xIang,devsq=0.;
	start=std::clock();iswrite=true;echeck1=false;echeck2=false;echeck3=false;
	kmin=4;kmax=10;//for tests with intensity
	testN=atoi(descr);ind=co;
int fmin=6850,fmax=9850,sep=(fmax-fmin)/(ind-1);//, ind=(fmax-fmin)/sep+1;//9950//11shots => sep=300

	if(sp==0){rhonor=1030713.25; heat=0.42107;  th=0.7328;fdiff=60;}//4.2802
	if(sp==1){rhonor= 767789.57; heat=0.42077;  th=1.9056;fdiff=60;}//2.4738//update Mar2012
	if(sp==2){rhonor= 817113.52; heat=0.37239;  th=2.0159;fdiff=60;}//3.0837
	if(sp==3){rhonor= 656228.94; heat=0.39849;  th=2.2074;fdiff=60;}//4.1912
	if(sp==4){rhonor= 397281.19; heat=0.41343;  th=2.1439;fdiff=60;}//3.0059

for(kk=kmin;kk<=kmax;kk++){xtotin[kk]=0.;xLPo[kk]=0.;xCP[kk]=0.;ytotin[kk]=0.;yLPo[kk]=0.;yCP[kk]=0.;};

start=std::clock();
for(fnum=fmin;fnum<=fmax;fnum+=sep){init(sp,fmin,fmax,sep);
		#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) shared(ittot)
		#include "intensity.cpp"
		for(kk=kmin;kk<=kmax;kk++){xtotin[kk]+=totin[kk]/ind;xLPo[kk]+=LPo[kk]/ind;xCP[kk]+=CP[kk]/ind;};};
		ans=( std::clock() - start) / (doub)CLOCKS_PER_SEC;//printf ("Time = %.2f s; \n", xans);

dheat=0.;drho=0.;dtheta=0.;

switch (testN){
case 1: snxy=75;break;
case 2: snxy=161;break;
case 3: ss=3e-3;break;
case 4:	ss=3e-2;break;
case 5:	fact=0.8;break;
case 6:	fact=1.2;break;
case 7:	echeck1=true;break;
case 8:	echeck2=true;break;
case 9: echeck3=true;break;
case 10: ind=11;sep=(fmax-fmin)/(ind-1); break;
case 11: ind=41;sep=(fmax-fmin)/(ind-1); break;
case 12: fdiff=90;break;
case 13: fdiff=40;break;
case 14: fdiff=0;break;
case 15: dheat=0.192;break;
case 16: dheat=0.096;break;
case 17: dheat=0.048;break;
case 18: dheat=0.024;break;
case 19: dheat=0.012;break;
case 20: dheat=0.006;break;
case 21: dheat=0.003;break;
case 22: dheat=0.0015;break;
case 23: drho=0.32;break;
case 24: drho=0.16;break;
case 25: drho=0.08;break;
case 26: drho=0.04;break;
case 27: drho=0.02;break;
case 28: drho=0.01;break;
case 29: drho=0.005;break;
case 30: drho=0.0025;break;
case 31: dtheta=0.08;break;
case 32: dtheta=0.04;break;
case 33: dtheta=0.02;break;
case 34: dtheta=0.01;break;
case 35: dtheta=0.005;break;
case 36: dtheta=0.0025;break;
case 37: dtheta=0.00125;break;
case 38: dtheta=0.000625;break;
	// 60->90 transition - to show it's below xisq<0.1
	// 40->60 transition - to show it's above xisq>0.1
	// 0->60  transition - to show that overall change is not extremely large
	// at the end perform all tests with ind=21; fdiff=60
}

heat*=(1+dheat);rhonor*=(1+drho);th*=(1+dtheta);

start=std::clock();
for(fnum=fmin;fnum<=fmax;fnum+=sep){init(sp,fmin,fmax,sep);
	#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) shared(ittot)
	#include "intensity.cpp"
for(kk=kmin;kk<=kmax;kk++){ytotin[kk]+=totin[kk]/ind;yLPo[kk]+=LPo[kk]/ind;yCP[kk]+=CP[kk]/ind;};
};
	xans=( std::clock() - start) / (doub)CLOCKS_PER_SEC;//printf ("Time = %.2f s; \n", xans);

//for(kk=kmin;kk<=kmax;kk++){xtotin[kk]/=ind;xLPo[kk]/=ind;xCP[kk]/=ind;ytotin[kk]/=ind;yLPo[kk]/=ind;yCP[kk]/=ind;};
devsq=0.;

for(kk=kmin;kk<=kmax;kk++)devsq+=(xtotin[kk]-ytotin[kk])*(xtotin[kk]-ytotin[kk])/dFnu[kk];
devsq+=(xCP[7]-yCP[7])*(xCP[7]-yCP[7])/dCP/dCP;
devsq+=(xCP[8]-yCP[8])*(xCP[8]-yCP[8])/dCP/dCP;
devsq+=(xLPo[4]-yLPo[4])*(xLPo[4]-yLPo[4])/dLP[0]/dLP[0];
devsq+=(xLPo[7]-yLPo[7])*(xLPo[7]-yLPo[7])/dLP[1]/dLP[1];
devsq+=(xLPo[8]-yLPo[8])*(xLPo[8]-yLPo[8])/dLP[2]/dLP[2];
devsq/=9.;

printf ("Spin %.2f, test N=%d, chisq= %.5f; first in %.2f s; second %.4f times faster \n", a, testN, devsq,ans,ans/xans);

stringstream ytr;ytr<<(int)100*a<<"in"<<ind<<"N"<<testN;string stra = ytr.str();
FILE * pFile; pFile = fopen ((dir+"testa"+stra+".dat").c_str(),"a");
fprintf(pFile,"%.2f %d %.5f\n",a, testN, devsq);
fclose(pFile);
;}
/*
changes in size for "current" list
old sizes => new sizes - 
0.053 chisq for spin 0
0.035 chisq for spin 0.5
0.015 chisq for spin 0.9

fact*=1.2
0.123 chisq for spin 0
0.091 chisq for spin 0.5
0.012 chisq for spin 0.9

fact*=0.8
0.014 chisq for spin 0
0.042 chisq for spin 0.5
0.021 chisq for spin 0.9
On other hand fact=0.6 - significantly underestimates the intensity

xaccur=3e-4 to xaccur = 3e-3=> 9e-4 xisq + 1.09 speed-up - spin 0.5
xaccurr=5e-4 to xaccurr=1e-2=> 9e-4 xisq + 1.50 speed-up - spin 0.5
xfact => huge changes at 87GHz - due to jet emission - we cut it
xss=1e-2                  => can't be improved
xsnxy=85  to   xsnxy=121  => 0.019 xisq +0.50 speed-up => 85 is good?
	xsnxy=111 to xsnxy=161=> 5e-3 xisq  +0.58 speed-up => 111 is the best
xstep=1e-2 to xstep=0.1   => 8e-4 xisq + 1.08 speed-up
xstep=-0.09               => can't be improved
xIint=1e-10 to xIint=2e-9 => 2e-4 xisq + 1.05 speed-up - spin 0.5;
xIang=1.                  => can't be improved


Current list:
xaccur=3e-3; //+1. accuracy of geodesics computation
xaccurr=1e-2;//+2. accuracy of radiative transfer integration
xfact=1.0;   //3. size of the integration region
xss=1e-2;    //+4. fractional distance from BH horizon to start integration
xsnxy=111;    //+5. grid N in the picture plane
xstep=0.1;   //+6. step control of geodesic computation
xsstep=-0.09;//+7. step control of radiative transfer computation
xIint=2e-9;  //+8. initial intensity along each ray
xIang=1.;    //+9. initial polarized phases along each ray

Hi-res list:
xaccur=3e-4; //+1. accuracy of geodesics computation
xaccurr=5e-4;//+2. accuracy of radiative transfer integration
xfact=1.0;   //3. size of the integration region
xss=3e-3;    //+4. fractional distance from BH horizon to start integration
xsnxy=285;    //+5. grid N in the picture plane
xstep=1e-2;   //+6. step control of geodesic computation
xsstep=-0.06;//+7. step control of radiative transfer computation
xIint=1e-10;  //+8. initial intensity along each ray
xIang=0.1;    //+9. initial polarized phases along each ray

Previous list:
xaccur=3e-4; //1. accuracy of geodesics computation
xaccurr=5e-4;//2. accuracy of radiative transfer integration
xfact=1.0;   //3. size of the integration region
xss=1e-2;    //4. fractional distance from BH horizon to start integration
xsnxy=85;    //5. grid N in the picture plane
xstep=1e-2;  //6. step control of geodesic computation
xsstep=-0.09;//7. step control of radiative transfer computation
xIint=1e-10; //8. initial intensity along each ray
xIang=1.;    //9. initial polarized phases along each ray


change 0.026 between "current" & "previous" list for spin 0   best
change 0.041 between "current" & "previous" list for spin 0.5 best
change 0.0028between "current" & "previous" list for spin 0.9 best

change 0.043 between "current" & "high-res" list for spin best 0
change 0.012 between "current" & "high-res" list for spin best 0.5
change 0.0025between "current" & "high-res" list for spin best 0.9

//----------------tests------------------------
	//if(echeck)rhopo+=0.2;//xisq=0.019
	//if(echeck)Upo-=0.1;//xisq=0.0076; Upo-=0.2 => xisq=0.04
	//if(echeck)u[1]*=temp;//xisq=0.00065;
	//if(echeck)u[2]*=temp;//xisq=0.00005;
	//if(echeck)u[3]*=temp;//xisq=0.00057;
	//if(echeck)Bi[1]*=temp;//xisq=0.0097;
	//if(echeck)Bi[2]*=temp;//xisq=0.0011;
	//if(echeck)Bi[3]*=temp;//xisq=0.0020;
//--------------end---tests--------------------
*/
/*
dheat=0.192; => xisq=6.67601;
dheat=0.096; => xisq=2.67034;
dheat=0.048; => xisq=0.69097;
dheat=0.024; => xisq=0.17334;
dheat=0.012; => xisq=0.04627;
dheat=0.006; => xisq=0.01321;
dheat=0.003; => xisq=0.00402;
dheat=0.0015; => xisq=0.00114;

drho=0.32; => xisq=;
drho=0.16; => xisq=1.49231;
drho=0.08; => xisq=0.56627;
drho=0.04; => xisq=0.17515;
drho=0.02; => xisq=0.04414;
drho=0.01; => xisq=0.01266;
drho=0.005; => xisq=0.00302;
drho=0.0025; => xisq=0.00154;

dtheta=0.08; => xisq=2.42186;
dtheta=0.04; => xisq=0.53909;
dtheta=0.02; => xisq=0.53535;
dtheta=0.01; => xisq=0.18599;
dtheta=0.005; => xisq=0.07712;
dtheta=0.0025; => xisq=0.20468;
dtheta=0.00125; => xisq=0.09617;
dtheta=0.000625; => xisq=0.01878;
*/