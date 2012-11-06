{doub xaccur,xaccurr,xfact,xss,xsnxy,xstep,xsstep,xIint,xIang,xans;
xaccur=3e-4; //+1. accuracy of geodesics computation
xaccurr=5e-4;//+2. accuracy of radiative transfer integration
xfact=1.0;   //3. size of the integration region
xss=3e-3;    //+4. fractional distance from BH horizon to start integration
xsnxy=285;    //+5. grid N in the picture plane
xstep=1e-2;   //+6. step control of geodesic computation
xsstep=-0.06;//+7. step control of radiative transfer computation
xIint=1e-10;  //+8. initial intensity along each ray
xIang=0.1;    //+9. initial polarized phases along each ray
accur=xaccur;accurr=xaccurr;fact=xfact;ss=xss;int snxy=xsnxy;doub step=xstep,sstep=xsstep,Iint=xIint,Iang=xIang,devsq=0.;
//! using h
	start=std::clock();iswrite=true;fact=1.0;co=0;
kmin=2;kmax=2;//for imaging

boosted=true;
if(boosted){sp=3;rhonor=527557.; heat=0.413558; th=2.234;boosted=true;} //0.9 spin best FnuLPCP fit // 4 2 1 4 __ 4 3 1 4 -- are command lines
if(!boosted){sp=3;rhonor=2475739.; heat=0.2015; th=2.653;boosted=false;}//0.9 spin best FnuLPCP fit - new unboosted// 4 2 0 4 __ 4 3 0 4 -- are command lines

//kmin=1;kmax=1;snxy=101;sp=1; rhonor=1476420.;heat=0.34531;th=2.29;boosted=false;
sp=3;rhonor=527557.; heat=0.413558; th=2.234;boosted=true;

printf("th=%f\n",th);tth=th;xth=th;sep=5;
	//for(fnum=-1-co;fnum<=-1-co;fnum+=sep){init(sp);//fnum<0 => averaged, fnum>0 => not averaged
	for(fnum=6950-co;fnum<=6950-co;fnum+=sep){init(sp);
	#pragma omp parallel for schedule(dynamic,1) shared(ittot)
	#include "imaging.cpp"
	ans=( std::clock() - start) / (doub)CLOCKS_PER_SEC;printf ("Time = %.2f s; th=%.3f; heat=%.3f\n", ans,th,heat);
	;}
}