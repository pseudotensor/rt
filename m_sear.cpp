{//search for minimum with a steepest descent method, less reliable than "m_space", but faster
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
	 devsq=0.;
doub resid[sflen][4], //residuals for \chi^2 minimization algorithm
	 Jac[sflen][3],   //correspondent Jacobian
	 bb[3],           //RHS vector for the next minimization iteration
	 matr[3][3],      //LHS matrix for the next minimization iteration
	 det;
doub ytotin[sflen][4],//for calculations of auxiliary spectrum
	 yLPo[sflen][4], 
	 yCP[sflen][4], 
	 yEVPA[sflen][4];
int snxy=xsnxy,       //global variable correspondent to xsnxy
    testN,            //particular deviation from given initial conditions - used to unrigorously avoid local minima in steepest descent
	ind,              //number of time frames for spectrum evaluation
	fmin,             //minimum ID of fluid simulation snapshot
	fmax,             //maximum ID of fluid simulation snapshot
	nP,               //maximum number of quantites to fit for
	nPeff,            //actual number of quantites we fit for
	oo;               //index of the set of parameters (rhonor,th,heat)

accur=xaccur;         //assigning values to global variables, which control radiative transfer
accurr=xaccurr;
fact=xfact;
ss=xss;

//selecting one model to test
sp=0;fmin=12424;fmax=22424;rhonor=16689876.83473; heat=0.53157;th=2.5056;thlimit=0.1;fdiff=0; isBcut=false;isBred=false;

ind=co;               //number of snapshots = 2nd command line argument
testN=atoi(descr);    //ID of initial deviation from global defined model = job array ID
nP=12;                //fit at most for 7 total fluxes, 3 LP fractions, 2 CP fractions
if(isLP87)            //either include or not include in the fit the LP fraction at 87GHz
	nPeff=nP; 
else 
	nPeff=nP-1;
sep=(fmax-fmin)/(ind-1);//compute difference of IDs of two consecutive snapshots
printf("testN=%d;  ind=%d\n",testN,ind);

dheat=0.;             //initialize local model (w/ deviation) = global model
drho=0.;
dtheta=0.;
switch (testN){       //define a particular local model based on value of testN
	case 1: break;											
	case 2: dheat=0.15;dheat=0.9;break;								
	case 3: dheat=-0.15;dheat=-0.5;break;								
	case 4:	drho=0.25;drho=0.9;break;						
	case 5:	drho=-0.25;drho=-0.6;break;
	case 6:	dtheta=0.1;dtheta=0.5;break;
	case 7:	dtheta=-0.1;dtheta=-0.5;break;
}
heat*=(1+dheat);      //compute model parameters for a local model
rhonor*=(1+drho);
th+=dtheta;

printf("in=%d shots; sp=%d; th=%.4f; heat=%.4f; rhonor=%.1f; fdiff=%d\n",co,sp,th,heat,rhonor,fdiff);
doub ddr=1.,          //differences of basic model parameters over 1 steepest descent iteration
	 ddh=0.02, 
	 dth=0.01;
niter=0;           
iswrite=false;        //do not write each computed model to file, since models are only slightly different by dheat, drho, dtheta

dheat=0.02;           //small deviations of basic model parameters for the purpose of computing a Jacobian
drho=0.02;            //these numbers are determined with consideration of convergence tests
dtheta=0.01;

while((fabs(ddh)>0.003)||(fabs(ddr)>0.01)||(fabs(dth)>0.005)){//convergence is stated if ALL model parameters change only slightly over 1 iteration
	niter++;          //iteration counter
	inp[0][0]=heat;   //0-th model + small deviations along "heat" axis, "rhonor" axis, and "theta" axis
	inp[0][1]=rhonor;
	inp[0][2]=th;
	
	inp[1][0]=heat*(1+dheat);
	inp[1][1]=rhonor;
	inp[1][2]=th;
	
	inp[2][0]=heat;
	inp[2][1]=rhonor*(1+drho);
	inp[2][2]=th;
	
	inp[3][0]=heat;
	inp[3][1]=rhonor;
	inp[3][2]=th*(1+dtheta);
	start=clock();
	
	//computing spectra for aforementioned 4 models
	for(oo=0;oo<4;oo++)//initialize these spectra with zeros
		for(kk=kmin;kk<=kmax;kk++){
			ytotin[kk][oo]=0.;
			yLPo[kk][oo]=0.;
			yCP[kk][oo]=0.;
			yEVPA[kk][oo]=0.;
		};
	for(fnum=fmin;fnum<=fmax;fnum+=sep){//loop over fluid simulation snapshots
		for(oo=0;oo<4;oo++){            //loop over 4 models
			heat=inp[oo][0];
			rhonor=inp[oo][1];
			th=inp[oo][2];
			init(sp,fmin,fmax,sep);     //compute spectrum for each fluid simulation snapshot for each model
			#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) shared(ittot)
			#include "intensity.cpp"
			for(kk=kmin;kk<=kmax;kk++){ //add spectra up to get average spectra
				ytotin[kk][oo]+=totin[kk]/ind;
				yLPo[kk][oo]+=LPo[kk]/ind;
				yCP[kk][oo]+=CP[kk]/ind;
				yEVPA[kk][oo]+=EVPA[kk]/ind;
			};
		};
	};
//end of integration block 
	
	for(oo=0;oo<4;oo++){                //compute normalized residuals with respect to the observed spectra
		for(il=4;il<=10;il++)
			resid[il-4][oo]=(ytotin[il][oo]-tofit[il][1])/dFnu[il];
		resid[7][oo]=(yLPo[4][oo]-tofit[4][2])/dLP[0];
		resid[8][oo]=(yLPo[7][oo]-tofit[7][2])/dLP[1];
		resid[9][oo]=(yLPo[8][oo]-tofit[8][2])/dLP[2];
		resid[10][oo]=(yCP[7][oo]-tofit[7][4])/dCP;
		resid[11][oo]=(yCP[8][oo]-tofit[8][4])/dCP;
	}
	{doub xisq=0.;
	for(il=0;il<nP;il++)                //calculate reduced \chi^2
		xisq+=resid[il][0]*resid[il][0]/(nPeff-3);
	printf("heat=%.4f, rhonor=%.1f, th=%.4f, xisq=%.3f\n",inp[0][0], inp[0][1], inp[0][2], xisq);

	stringstream sss;                   //write into "xisqa" file 0-th model parameters and correspondent reduced \chi^2
	sss<<(int)100*a<<"in"<<co<<"N"<<testN<<"fdiff"<<fdiff;
	string outstr=dir+"xisqa"+sss.str(); 
	outstr+=".dat";
	FILE * xFile; 
	xFile = fopen (outstr.c_str(),"a");
	fprintf(xFile,"%.2f %.4f %.4f %.5f %.5f %.3e\n",a,xisq,inp[0][2],inp[0][0],inp[0][1],rate*year/Msun);
	fclose(xFile);
	};

	ans=(clock() - start ) / (doub)CLOCKS_PER_SEC;//timing
	printf ("Time = %.2f s; finished %d iterations; th=%.3f; heat=%.3f\n", ans,niter,th,heat);
	int ij, ii;
	for(il=0;il<nP;il++){                         //compute a Jacobian based on differences in residuals
		Jac[il][0]=(resid[il][1]-resid[il][0])/dheat;
		Jac[il][1]=(resid[il][2]-resid[il][0])/drho;
		Jac[il][2]=(resid[il][3]-resid[il][0])/dtheta;
	};
	for(ij=0;ij<3;ij++)                          //compute LHS matrix as Jacobian^2 (w/ transposition)
		for(ii=0;ii<3;ii++){
			matr[ij][ii]=0.;
			for(il=0;il<nP;il++)
				matr[ij][ii]+=Jac[il][ij]*Jac[il][ii];
		};

	for(il=0;il<3;il++)
		matr[il][il]*=1.1;                       //multiply diagonal elements of LHS matrix by some number following Levenberg-Marquardt algorithm - may want to play with that!

	for(ii=0;ii<3;ii++){                         //calculate RHS vector of steepest descent method
		bb[ii]=0.;
		for(il=0;il<nP;il++)
			bb[ii]+=Jac[il][ii]*resid[il][0];
	};                                           //compute a determinant and solve a linear system of equations for parameters deviations for the next iteration
	det=-(matr[0][2]*matr[1][1]*matr[2][0]) + matr[0][1]*matr[1][2]*matr[2][0] + matr[0][2]*matr[1][0]*matr[2][1] - matr[0][0]*matr[1][2]*matr[2][1] - matr[0][1]*matr[1][0]*matr[2][2] + matr[0][0]*matr[1][1]*matr[2][2];
	ddh=(bb[2]*matr[0][2]*matr[1][1] - bb[2]*matr[0][1]*matr[1][2] - bb[1]*matr[0][2]*matr[2][1] + bb[0]*matr[1][2]*matr[2][1] + (bb[1]*matr[0][1] - bb[0]*matr[1][1])*matr[2][2])/det;
	ddr=(bb[2]*(-(matr[0][2]*matr[1][0]) + matr[0][0]*matr[1][2]) + bb[1]*(matr[0][2]*matr[2][0] - matr[0][0]*matr[2][2]) + bb[0]*(-(matr[1][2]*matr[2][0]) + matr[1][0]*matr[2][2]))/det;
	dth=(bb[2]*matr[0][1]*matr[1][0] - bb[2]*matr[0][0]*matr[1][1] - bb[1]*matr[0][1]*matr[2][0] + bb[0]*matr[1][1]*matr[2][0] + (bb[1]*matr[0][0] - bb[0]*matr[1][0])*matr[2][1])/det;
	printf("Iterated ddr=%.4f, ddh=%.4f, dth=%.4f\n",ddr,ddh,dth);
//done till here
	while((fabs(ddh)>0.10) || (fabs(ddr)>0.2) || (fabs(dth)>0.05)){//while changes over 1 iteration are too large, then halve these changes
		ddh/=2.;
		ddr/=2.;
		dth/=2.;
	};
	doub xf=1.0;                                //update by an "xf" fraction of a vector
	heat=inp[0][0]*(1+xf*ddh);                  //update quantities for the next iteration
	rhonor=inp[0][1]*(1+xf*ddr);
	th=inp[0][2]*(1+xf*dth);
	if (niter>20){                              //limit to 20 iterations
		break;
	};
};


//compute spectrum for the final best-fitting solution
stringstream sss;sss<<(int)100*a;
//beginning of integration block
for(kk=kmin;kk<=kmax;kk++){
	xtotin[kk]=0.;
	xLPo[kk]=0.;
	xCP[kk]=0.;
	xEVPA[kk]=0.;
};
for(fnum=fmin;fnum<=fmax;fnum+=sep){            //cycle over all fluid simulation snapshots
	init(sp,fmin,fmax,sep);
	#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) shared(ittot)
	#include "intensity.cpp"
	for(kk=kmin;kk<=kmax;kk++){                 //compute the spectrum
		xtotin[kk]+=totin[kk]/ind;
		xLPo[kk]+=LPo[kk]/ind;
		xCP[kk]+=CP[kk]/ind;
		xEVPA[kk]+=EVPA[kk]/ind;
	};
}
//end of integration block 	
for(il=4;il<=10;il++)                           //compute residuals
	resid[il-4][0]=(xtotin[il]-tofit[il][1])/dFnu[il];
resid[7][oo]=(xLPo[4]-tofit[4][2])/dLP[0];
resid[8][oo]=(xLPo[7]-tofit[7][2])/dLP[1];
resid[9][oo]=(xLPo[8]-tofit[8][2])/dLP[2];
resid[10][oo]=(xCP[7]-tofit[7][4])/dCP;
resid[11][oo]=(xCP[8]-tofit[8][4])/dCP;
doub xisq=0.;
for(il=0;il<nP;il++)                            //compute reduced \chi^2
	xisq+=resid[il][0]*resid[il][0]/(nPeff-3);

string stra = sss.str();                        //write best-fitting model into "bestfita" file
FILE * pFile;
pFile = fopen ((dir+"bestfita"+stra+".dat").c_str(),"a");
for(kk=kmin;kk<=kmax;kk++){
	printf("avg at f=%.1f; I=%.3fJy LP=%.2f%% CP=%.3f%% EVPA=%.2fdeg\n",sftab[kk][0], xtotin[kk],xLPo[kk], xCP[kk],xEVPA[kk]);
	fprintf(pFile,"%.1f %.3f %.3f %.3f %.2f %.5f %.4f %.4f %.4f\n",sftab[kk][0], xtotin[kk],xLPo[kk], xCP[kk],xEVPA[kk],th,heat,rhonor,Bpo);
};
fclose(pFile);
}