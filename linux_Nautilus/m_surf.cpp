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

//for(kk=kmin;kk<=kmax;kk++){xtotin[kk]=0.;xLPo[kk]=0.;xCP[kk]=0.;xEVPA[kk]=0.;};
//0.5 spin best FnuLPCP fit // 2 2 2.0 3 __ 2 3 2.0 3 -- are command lines; PBS_ARRAYID=0..10
	if(sp==0){rhonor=1030713.25; heat=0.42107;  th=0.7328;}//4.2802
	if(sp==1){rhonor= 988988.82; heat=0.37012;  th=1.8407;}//2.9363
	if(sp==2){rhonor= 817113.52; heat=0.37239;  th=2.0159;}//3.0837
	if(sp==3){rhonor= 656228.94; heat=0.39849;  th=2.2074;}//4.1912
	if(sp==4){rhonor= 397281.19; heat=0.41343;  th=2.1439;}//3.0059

int nnx=40;
iswrite=true;oo=0;
//doub chmax=0.045,crmax=0.09,cthmax=0.15;//for spin 0.5
doub chmax=0.045,crmax=0.09,cthmax=0.1;//for spin 0.9
kmin=4;kmax=10;cheat=heat; crhonor=rhonor;doub ch=0.,cr=0.,delta=1e-10, xth=th;

int ind=atoi(descr),fmin=6950,fmax=9950,sep=300;//ind=0..10
fnum=fmin+ind*sep;

if(co==1){// accretion rate & theta - not used for plots
	fif="cr_th";
for(cr=-crmax;cr<=crmax;cr+=2.*crmax/(nnx-1)-delta)
for(tth=xth-cthmax;tth<=xth+cthmax;tth+=2.*cthmax/(nnx-1)-delta){
heat=cheat*(1+ch);rhonor=crhonor*(1+cr);th=tth;

{init(sp,fmin,fmax,sep);start=std::clock();
	#pragma omp parallel for schedule(dynamic,1) shared(ittot)
	#include "intensity.cpp"
}

;};}

if(co==2){// accretion rate & heating constant - used for first plot
	fif="ch_cr";
for(ch=-chmax;ch<=chmax;ch+=2.*chmax/(nnx-1)-delta)
for(cr=-crmax;cr<=crmax;cr+=2.*crmax/(nnx-1)-delta){
//for(th=tth-cthmax;th<=tth+cthmax;th+=2.*cthmax/(nnx-1)-delta){
heat=cheat*(1+ch);rhonor=crhonor*(1+cr);

{init(sp,fmin,fmax,sep);start=std::clock();
	#pragma omp parallel for schedule(dynamic,1) shared(ittot)
	#include "intensity.cpp"
}

;};}

if(co==3){// heating constant & theta - used for second plot
	fif="ch_th";
for(ch=-chmax;ch<=chmax;ch+=2.*chmax/(nnx-1)-delta)
//for(cr=-crmax;cr<=crmax;cr+=2.*crmax/(nnx-1)-delta)
for(tth=xth-cthmax;tth<=xth+cthmax;tth+=2.*cthmax/(nnx-1)-delta){
heat=cheat*(1+ch);rhonor=crhonor*(1+cr);th=tth;

{init(sp,fmin,fmax,sep);start=std::clock();
	#pragma omp parallel for schedule(dynamic,1) shared(ittot)
	#include "intensity.cpp"
}


;}
;}
;}