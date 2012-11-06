{doub xaccur,xaccurr,xfact,xss,xsnxy,xstep,xsstep,xIint,xIang,xans;
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

	start=std::clock();sep=50;boosted=true;iswrite=true;fact=1.0;co=0;

kmin=4;kmax=10;//for tests with intensity
    heat=0.33150;rhonor=575101.7855;th=PI/2.*(2*cas-1)/thn;
sp=1; rhonor=640000.;heat=0.39; th=2.217; boosted=false;

sp=1; rhonor=300090.;heat=0.443782; th=1.86109; boosted=true;//for tests

//sp=0; rhonor=792833.; heat=0.439559;  th=2.10488;boosted=true;
//sp=1;rhonor=800628.; heat=0.434947;  th=2.00032;boosted=true;//new
sp=2; rhonor=462221.; heat=0.470756;  th=2.04975;boosted=true;
//sp=3;rhonor=760871.; heat=0.379966;  th=2.18908;boosted=false;
//sp=3;rhonor=527557.; heat=0.413558;  th=2.234;boosted=true;//larger coeff at rVc => larger EVPA; B[r]<0 along the ray approaching observer, rVc<0
//sp=4; rhonor=304073.; heat=0.423505;  th=2.20845;boosted=true;
//sp=4; rhonor=6e5;th=2.14674;// tests with Lei

	printf("th=%f\n",th);tth=th;xth=th;
	for(fnum=-1-co;fnum<=-1-co;fnum+=sep){init(sp);//fnum<0 => averaged, fnum>0 => not averaged
	//for(fnum=6950-co;fnum<=6950-co;fnum+=sep){init(sp);
	#pragma omp parallel for schedule(dynamic,1) shared(ittot)
	#include "intensity.cpp"
	inited=false;init(sp);
	ans=( std::clock() - start) / (doub)CLOCKS_PER_SEC;printf ("Time = %.2f s; th=%.3f; heat=%.3f\n", ans,th,heat);
	;}
}