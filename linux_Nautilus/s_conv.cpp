{doub xaccur,xaccurr,xfact,xss,xsnxy,xstep,xsstep,xIint,xIang,xans;
xaccur=3e-3; //+1. accuracy of geodesics computation
xaccurr=1e-2;//+2. accuracy of radiative transfer integration
xfact=1.0;   //3. size of the integration region
xss=1e-2;    //+4. fractional distance from BH horizon to start integration
xsnxy=111;    //+5. grid N in the picture plane
xstep=0.1;   //+6. step control of geodesic computation
xsstep=-0.09;//+7. step control of radiative transfer computation
xIint=2e-9;  //+8. initial intensity along each ray
xIang=1.;    //+9. initial polarized phases along each ray
accur=xaccur;accurr=xaccurr;fact=xfact;ss=xss;int snxy=xsnxy;doub step=xstep,sstep=xsstep,Iint=xIint,Iang=xIang,devsq=0.;
//echeck=false;
	start=std::clock();sep=50;iswrite=true;co=0;th=PI/2.*(2*cas-1)/thn;
	kmin=4;kmax=10;//for tests with intensity

//sp=0; rhonor=792833.; heat=0.439559;  th=2.10488;boosted=true;
sp=1;rhonor=800628.; heat=0.434947;  th=2.00032;boosted=true;//new
//sp=2; rhonor=462221.; heat=0.470756;  th=2.04975;boosted=true;
//sp=3;rhonor=760871.; heat=0.379966;  th=2.18908;boosted=false;
//sp=3;rhonor=527557.; heat=0.413558;  th=2.234;boosted=true;
//sp=4; rhonor=304073.; heat=0.423505;  th=2.20845;boosted=true;



	fnum=-1;//fnum=6950;
	init(sp);{start=std::clock();
		#pragma omp parallel for schedule(dynamic,1) shared(ittot)
		#include "intensity.cpp"
	};xans=( std::clock() - start) / (doub)CLOCKS_PER_SEC;//printf ("Time = %.2f s; \n", xans);

//snxy=75;//test 1
//snxy=161;//test 2
//ss=3e-2;//test 3
//fact=0.8;//test 4
//test 5
//echeck=true;
//fact=1.2;
echeck3=true;

init(sp);

	for(kk=kmin;kk<=kmax;kk++){xtotin[kk]=totin[kk];xLPo[kk]=LPo[kk];xCP[kk]=CP[kk];}{start=std::clock();
		#pragma omp parallel for schedule(dynamic,1) shared(ittot)
		#include "intensity.cpp"
	};devsq=0.;

for(kk=kmin;kk<=kmax;kk++)devsq+=(xtotin[kk]-totin[kk])*(xtotin[kk]-totin[kk])/dFnu[kk];
devsq+=(xCP[7]-CP[7])*(xCP[7]-CP[7])/dCP/dCP;
devsq+=(xCP[8]-CP[8])*(xCP[8]-CP[8])/dCP/dCP;
devsq+=(log(xLPo[4])-log(LPo[4]))*(log(xLPo[4])-log(LPo[4]))/dLP[0]/dLP[0];
devsq+=(log(xLPo[7])-log(LPo[7]))*(log(xLPo[7])-log(LPo[7]))/dLP[1]/dLP[1];
devsq+=(log(xLPo[8])-log(LPo[8]))*(log(xLPo[8])-log(LPo[8]))/dLP[2]/dLP[2];
devsq/=9.;

ans=( std::clock() - start) / (doub)CLOCKS_PER_SEC;//printf ("Time = %.2f s; \n", ans);
printf ("chisq= %.5f; first in %.2f s; second %.4f times faster \n", devsq,xans,xans/ans);
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