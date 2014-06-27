{//surf region close to best fit
doub xaccur=3e-3,  //1. absolute accuracy of geodesics computation
	 xaccurr=1e-2, //2. absolute accuracy of radiative transfer integration
     xfact=1.0,    //3. relative size of the integration region
	 xss=1e-2,     //4. fractional distance from BH horizon to the sphere, where geodesic integration stops
	 xsnxy=111,    //5. number of points N along one side in the picture plane for N x N intensity grid
	 xstep=0.1,    //6. step size in geodesic computation
	 xsstep=-0.09, //7. step size in radiative transfer computation
	 xIint=2e-9,   //8. initial intensity along each ray for radiative transfer
	 xIang=1.;     //9. initial polarized phases along each ray for radiative transfer
doub step=xstep,   //local variables, which control radiative transfer
	 sstep=xsstep,
	 Iint=xIint,
	 Iang=xIang;
int snxy=xsnxy;    //global variable correspondent to xsnxy

accur=xaccur;      //assigning values to global variables, which control radiative transfer
accurr=xaccurr;
fact=xfact;
ss=xss;

//must specify models for each spin - choose your models!
if(sp==0){rhonor=1030713.25; heat=0.42107;  th=0.7328;}
if(sp==1){rhonor= 988988.82; heat=0.37012;  th=1.8407;}
if(sp==2){rhonor= 817113.52; heat=0.37239;  th=2.0159;}
if(sp==3){rhonor= 656228.94; heat=0.39849;  th=2.2074;}
if(sp==4){rhonor= 397281.19; heat=0.41343;  th=2.1439;}

int nnx=40;         //number of steps along each variable for models slightly different from the best-fitting
doub cheat=heat,    //auxiliary heating constant, density normalization, and viewing angle
	 crhonor=rhonor,
	 ctheta=th;
iswrite=true;       //write to file
doub chmax=0.045,   //maximum deviations along each variable from the best-fitting model
	 crmax=0.09,
	 cthmax=0.1;
kmin=4; kmax=10;    //standard subset of frequencies for spectrum evaluation
doub ch=0.,         //step variables
	 cr=0.,
	 cth=0.,
	 delta=1e-10,   //small delta ensures last value of step variable is still within the range
	 xth=th;

int ind=atoi(descr),//job array ID = number of fluid simulation snapshots
	fmin=6950,      //specify the minimum and maximum fluid simulation snapshot ID
	fmax=9950,
	sep=(fmax-fmin)/(ind-1);//calculate ID difference between consecutive snapshots

//three man cycles - surfing three planes of two pairs of model parameters
if(co==1){// accretion rate & viewing angle theta
	fif="cr_th";
	for(cr=-crmax;cr<=crmax;cr+=2.*crmax/(nnx-1)-delta)
		for(cth=-cthmax;cth<=cthmax;cth+=2.*cthmax/(nnx-1)-delta){
			heat=cheat*(1+ch);
			rhonor=crhonor*(1+cr);
			th=ctheta+cth;
			for(fnum=fmin;fnum<=fmax;fnum+=sep){             //cycle over fluid simulation snapshots - calculating individual spectra. Average spectra & \chi^2 are computed by separate Mathematica scripts
				init(sp,fmin,fmax,sep);
				start=clock();
				#pragma omp parallel for schedule(dynamic,1) shared(ittot)
				#include "intensity.cpp"
			};
		};
};

if(co==2){// accretion rate & heating constant
	fif="ch_cr";
	for(ch=-chmax;ch<=chmax;ch+=2.*chmax/(nnx-1)-delta)
		for(cr=-crmax;cr<=crmax;cr+=2.*crmax/(nnx-1)-delta){
			heat=cheat*(1+ch);
			rhonor=crhonor*(1+cr);
			th=ctheta+cth;
			for(fnum=fmin;fnum<=fmax;fnum+=sep){             //cycle over fluid simulation snapshots
				init(sp,fmin,fmax,sep);
				start=clock();
				#pragma omp parallel for schedule(dynamic,1) shared(ittot)
				#include "intensity.cpp"
			};
		};
};

if(co==3){// heating constant & viewing angle theta
	fif="ch_th";
	for(ch=-chmax;ch<=chmax;ch+=2.*chmax/(nnx-1)-delta)
		for(cth=-cthmax;cth<=cthmax;th+=2.*cthmax/(nnx-1)-delta){
			heat=cheat*(1+ch);
			rhonor=crhonor*(1+cr);
			th=ctheta+cth;
			for(fnum=fmin;fnum<=fmax;fnum+=sep){             //cycle over fluid simulation snapshots
				init(sp,fmin,fmax,sep);
				start=clock();
				#pragma omp parallel for schedule(dynamic,1) shared(ittot)
				#include "intensity.cpp"
			};
		};
};
};