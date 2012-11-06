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
//sp=3;heat=0.345315;rhonor=643568.;th=2.14674;//0.9 spin best FnuCPLP fit
		//sp=1;heat=0.359703;rhonor=466429.;th=1.93731;//0.5 spin best FnuCPLP fit
	//sp=3;heat=0.33150;rhonor=818349.3297; th=2.0420; boosted=false;//0.9 spin best FnuLPCP fit - new
if(boosted && (sp==1)){rhonor=297000.;heat=0.448; th=1.893;} //0.5 spin best FnuLPCP fit // 2 2 1 3 __ 2 3 1 3 -- are command lines
if((!boosted) && (sp==1)){rhonor=640000.;heat=0.39; th=2.217;}//0.9 spin best FnuLPCP fit - new boosted// 2 2 0 3 __ 2 3 0 3 -- are command lines

if(boosted && (sp==3)){rhonor=527557.; heat=0.413558;  th=2.234;} 
if((!boosted) && (sp==3)){rhonor=760871.; heat=0.379966;  th=2.18908;}

int nnx=9;
iswrite=true;oo=0;//doub chmax=0.045,crmax=0.09,cthmax=0.15;//for spin 0.5
doub chmax=0.045,crmax=0.09,cthmax=0.1;//for spin 0.9
kmin=4;kmax=10;cheat=heat; crhonor=rhonor;doub ch=0.,cr=0.,delta=1e-10, xth=th;

if(co==1){// accretion rate & theta - not used for plots
	fif="cr_th";
for(cr=-crmax;cr<=crmax;cr+=2.*crmax/(nnx-1)-delta)
for(tth=xth-cthmax;tth<=xth+cthmax;tth+=2.*cthmax/(nnx-1)-delta){
heat=cheat*(1+ch);rhonor=crhonor*(1+cr);th=tth;
for(fnum=-1;fnum<=-1;fnum+=sep){init(sp);//4379->5699
start=std::clock();
#pragma omp parallel for schedule(dynamic,1) shared(ittot)
#include "intensity.cpp"
;};};}

if(co==2){// accretion rate & heating constant - used for first plot
	fif="ch_cr";
for(ch=-chmax;ch<=chmax;ch+=2.*chmax/(nnx-1)-delta)
for(cr=-crmax;cr<=crmax;cr+=2.*crmax/(nnx-1)-delta){
//for(th=tth-cthmax;th<=tth+cthmax;th+=2.*cthmax/(nnx-1)-delta){
heat=cheat*(1+ch);rhonor=crhonor*(1+cr);
for(fnum=-1;fnum<=-1;fnum+=sep){init(sp);//4379->5699
start=std::clock();
#pragma omp parallel for schedule(dynamic,1) shared(ittot) firstprivate(a)
#include "intensity.cpp"
;};};}

if(co==3){// heating constant & theta - used for second plot
	fif="ch_th";
for(ch=-chmax;ch<=chmax;ch+=2.*chmax/(nnx-1)-delta)
//for(cr=-crmax;cr<=crmax;cr+=2.*crmax/(nnx-1)-delta)
for(tth=xth-cthmax;tth<=xth+cthmax;tth+=2.*cthmax/(nnx-1)-delta){
heat=cheat*(1+ch);rhonor=crhonor*(1+cr);th=tth;
for(fnum=-1;fnum<=-1;fnum+=sep){init(sp);//4379->5699
start=std::clock();
#pragma omp parallel for schedule(dynamic,1) shared(ittot) firstprivate(a)
#include "intensity.cpp"
;}
;}
;}
;}