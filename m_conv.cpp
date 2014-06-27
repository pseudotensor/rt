{//convergence studies
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
	 Iang=xIang,
	 chisq=0.,     //\chi^2/dof
	 xans;         //calculation time for modified model (to compare with calculation time of original model)
doub ztotin[sflen],//spectrum for the modified model
	 zLPo[sflen], 
	 zCP[sflen];
int snxy=xsnxy,    //global variable correspondent to xsnxy
    testN,         //test number
	ind;           //number of time frames for spectrum evaluation

accur=xaccur;      //assigning values to global variables, which control radiative transfer
accurr=xaccurr;
fact=xfact;
ss=xss;
start=clock();
iswrite=true;      //write output to file
echeck1=false;     //switch off checks for the unmodified model
echeck2=false;
echeck3=false;

kmin=4;kmax=10;    //standard set of frequencies for tests
testN=atoi(descr); //job array ID = test number
ind=co;
int fmin=6850,     //default range of fluid simulation IDs
	fmax=9850,
	sep=(fmax-fmin)/(ind-1);

//must specify models for each spin - choose your models!
if(sp==0){rhonor=1030713.25; heat=0.42107;  th=0.7328;fdiff=60;}
if(sp==1){rhonor= 767789.57; heat=0.42077;  th=1.9056;fdiff=60;}
if(sp==2){rhonor= 817113.52; heat=0.37239;  th=2.0159;fdiff=60;}
if(sp==3){rhonor= 656228.94; heat=0.39849;  th=2.2074;fdiff=60;}
if(sp==4){rhonor= 397281.19; heat=0.41343;  th=2.1439;fdiff=60;}

for(kk=kmin;kk<=kmax;kk++){//initializing spectra with zeros
	xtotin[kk]=0.;
	xLPo[kk]=0.;
	xCP[kk]=0.;
	ztotin[kk]=0.;
	zLPo[kk]=0.;
	zCP[kk]=0.;
};

start=clock();
for(fnum=fmin;fnum<=fmax;fnum+=sep){             //cycle over fluid simulation snapshots - computing average spectrum
	init(sp,fmin,fmax,sep);
	#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) shared(ittot)
	#include "intensity.cpp"
	for(kk=kmin;kk<=kmax;kk++){
		xtotin[kk]+=totin[kk]/ind;
		xLPo[kk]+=LPo[kk]/ind;
		xCP[kk]+=CP[kk]/ind;
	};
};
ans=(clock() - start) / (doub)CLOCKS_PER_SEC;   //timing of unmodified model calculation

dheat=0.;                                       //set default modification for the model at zeros
drho=0.;
dtheta=0.;

switch (testN){                                 //specifying modified model
	case 1: snxy=75;break;                      //smaller and larger resolution in the picture plane
	case 2: snxy=161;break;
	case 3: ss=3e-3;break;                      //smaller and larger fractional distance from BH horizon to the sphere, where geodesic integration stops
	case 4:	ss=3e-2;break;
	case 5:	fact=0.8;break;                     //smaller and larger relative size of the integration region
	case 6:	fact=1.2;break;
	case 7:	echeck1=true;break;                 //changes to the power-laws of radial extension of quantities
	case 8:	echeck2=true;break;
	case 9: echeck3=true;break;
	case 10: ind=11;sep=(fmax-fmin)/(ind-1); break;//smaller and larger number of fluid simulation snapshots over which to compute average spectrum
	case 11: ind=41;sep=(fmax-fmin)/(ind-1); break;
	case 12: fdiff=90;break;                    //smaller and larger number of fluid simulation snapshots, over which the geodesics are propagated as the simulation is evolved
	case 13: fdiff=40;break;
	case 14: fdiff=0;break;                     //fast light approximation
	case 15: dheat=0.192;break;                 //range of relative changes to heating constant
	case 16: dheat=0.096;break;
	case 17: dheat=0.048;break;
	case 18: dheat=0.024;break;
	case 19: dheat=0.012;break;
	case 20: dheat=0.006;break;
	case 21: dheat=0.003;break;
	case 22: dheat=0.0015;break;
	case 23: drho=0.32;break;                   //range of relative changes to heating density
	case 24: drho=0.16;break;
	case 25: drho=0.08;break;
	case 26: drho=0.04;break;
	case 27: drho=0.02;break;
	case 28: drho=0.01;break;
	case 29: drho=0.005;break;
	case 30: drho=0.0025;break;
	case 31: dtheta=0.08;break;                 //range of relative (not absolute) changes to viewing angle
	case 32: dtheta=0.04;break;
	case 33: dtheta=0.02;break;
	case 34: dtheta=0.01;break;
	case 35: dtheta=0.005;break;
	case 36: dtheta=0.0025;break;
	case 37: dtheta=0.00125;break;
	case 38: dtheta=0.000625;break;
}

heat*=(1+dheat);                               //modifying the basic model parameters
rhonor*=(1+drho);
th*=(1+dtheta);

start=clock();
for(fnum=fmin;fnum<=fmax;fnum+=sep){           //cycle over fluid simulation snapshots - computing average spectrum for modified model
	init(sp,fmin,fmax,sep);
	#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) shared(ittot)
	#include "intensity.cpp"
	for(kk=kmin;kk<=kmax;kk++){
		ztotin[kk]+=totin[kk]/ind;
		zLPo[kk]+=LPo[kk]/ind;
		zCP[kk]+=CP[kk]/ind;
	};
};
xans=(clock() - start) / (doub)CLOCKS_PER_SEC;//timing the modified model calculation


for(kk=kmin;kk<=kmax;kk++)                    //computing effective \chi^2 between modified and unmodified model
	chisq+=(xtotin[kk]-ztotin[kk])*(xtotin[kk]-ztotin[kk])/dFnu[kk];
chisq+=(xCP[7]-zCP[7])*(xCP[7]-zCP[7])/dCP/dCP;
chisq+=(xCP[8]-zCP[8])*(xCP[8]-zCP[8])/dCP/dCP;
chisq+=(xLPo[4]-zLPo[4])*(xLPo[4]-zLPo[4])/dLP[0]/dLP[0];
chisq+=(xLPo[7]-zLPo[7])*(xLPo[7]-zLPo[7])/dLP[1]/dLP[1];
chisq+=(xLPo[8]-zLPo[8])*(xLPo[8]-zLPo[8])/dLP[2]/dLP[2];
chisq/=9.;                                    //switching to reduced \chi^2
printf ("Spin %.2f, test N=%d, chisq= %.5f; first in %.2f s; second %.4f times faster \n", a, testN, chisq,ans,ans/xans); //printing results on screen including comparison of runtimes

stringstream ytr;                             //printing the resuls of tests into "testa" file
ytr<<(int)100*a<<"in"<<ind<<"N"<<testN;
string stra = ytr.str();
FILE * pFile; 
pFile = fopen ((dir+"testa"+stra+".dat").c_str(),"a");
fprintf(pFile,"%.2f %d %.5f\n",a, testN, chisq);
fclose(pFile);
}