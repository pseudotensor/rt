{//surf the entire parameter space searching for the best fit to polarized spectrum
  //doub xaccur=3e-3,  //1. absolute accuracy of geodesics computation
 //RG: TRY TO AVOID: "Radial coordinate on a geodesic rr=0.822507 is out of bounds"

  doub xaccur=3e-3,  //1. absolute accuracy of geodesics computation
       xaccurr=1e-2, //2. absolute accuracy of radiative transfer integration
     xfact=1.0,    //3. relative size of the integration region
	 xss=1e-2,     //4. fractional distance from BH horizon to the sphere, where geodesic integration stops
	 xsnxy=101,    //5. number of points N along one side in the picture plane for N x N intensity grid
	 xstep=0.1,    //6. step size in geodesic computation
	 xsstep=-0.09, //7. step size in radiative transfer computation
	 xIint=2e-9,   //8. initial intensity along each ray for radiative transfer
	 xIang=1.;     //9. initial polarized phases along each ray for radiative transfer
doub step=xstep,   //local variables, which control radiative transfer
	 sstep=xsstep, 
  Iint=xIint,
  Iang=xIang;

 int snxy=xsnxy,    //global variable correspondent to xsnxy
 	ind=co,        //number of snapshots over which to compute mean intensity = 2-nd command line argument
 	oo,            //index of the set of parameters (rhonor,th,heat)
    //  fmin=2150,     //minimum ID of fluid simulation snapshot
 	// fmax=2150,     //maximum ID of fluid simulation snapshot
     fmin=5500,     //minimum ID of fluid simulation snapshot
 	fmax=5520,     //maximum ID of fluid simulation snapshot
 	sep=(fmax-fmin)/(ind-1);//ID difference between consecutive considered fluid simulation snapshots
 // sep =1;

doub resid[sflen][4], //residuals for \chi^2 minimization algorithm
	 Jac[sflen][3],   //correspondent Jacobian
	 bb[3],           //RHS vector for the next minimization iteration
	 matr[3][3];      //LHS matrix for the next minimization iteration

// #include "lightup_jet.cpp"

accur=xaccur;      //assigning values to global variables, which control radiative transfer
accurr=xaccurr;
fact=xfact;
ss=xss;

cas=atoi(descr);      //batch job ID = environment variable

//RG: DEFAULT: ASTRORAYv1
//th=PI/2.*(2*cas-1)/thn;//define theta angle out of a discrete set
//th=PI/2.-0.4;

start=clock();

rhonor=3e5;           //initial density for iterations
drho=0.02;            //density relative offset to compute Jacobian

// RG: can't reinitialize...
// #include "lightup_jet.cpp"
// fmin=5500;fmax=5520;sep=21;kmin=0;kmax=10;sp=0;rhonor=80000.; heat=0.0601; th=1.4; dphi=4.*PI/3.;thlimit=0.05;fdiff=20;isBcut=false;isBred=true;magn_cap=4;Te_jet_par=10.; include_jet=1;

// case 72000: 
// fmin=5500;fmax=5520;sep=21;kmin=4;kmax=10;sp=0;rhonor=371131.52176; heat=0.24979; th=0.6457718232379019; dphi=4.*PI/3.;thlimit=0.05;fdiff=20;isBcut=false;isBred=true;magn_cap=4;Te_jet_par=10.; include_jet=1;
// case 72001: 
// fmin=5500;fmax=5520;kmin=4;kmax=10;sp=0;rhonor=500000.; heat=0.24979; th=1.4; dphi=4.*PI/3.;thlimit=0.05;fdiff=20;isBcut=false;isBred=true;magn_cap=4;Te_jet_par=10.; include_jet=1;

// case 70006: 
fmin=5500;fmax=5600;sep=5;kmin=4;kmax=10;sp=0;rhonor=120304.04439; heat=0.15599; th=0.646; dphi=4.*PI/3.;thlimit=0.05;fdiff=20;isBcut=false;isBred=true;magn_cap=4;Te_jet_par=35.; include_jet=0; 


sep=(fmax-fmin)/(ind-1);//ID difference between consecutive considered fluid simulation snapshots

printf(YELLOW"[m_space.cpp]: "RESET"in=%d shots; sp=%d; angle th=%f\n",co,sp,th);

 doub heat_min=0.4;
// RG shorten it?
// for(heat=0.4;heat>0.39;heat*=0.96){             //choose heat=0.4
// for(heat=0.75;heat>0.15;heat*=0.96){               //choose a set of heat parameter
for(heat=heat_min;heat>0.1;heat*=0.96){               //choose a set of heat parameter

	doub ddr=1.;                                   //relative adjustment of density for next iteration
	niter=0;                                       //number of iterations
	iswrite=false;
	while(fabs(ddr)>0.01){                         //do iterations till convergence is reached
		niter++;

		 //initialize 2 states with same heat and slightly different rhonor
		inp[0][0]=heat;                          
		inp[0][1]=rhonor;
		inp[1][0]=heat;
		inp[1][1]=rhonor*(1+drho);

		//compute spectra for both states
		for(oo=0;oo<2;oo++){                      
			heat=inp[oo][0];
			rhonor=inp[oo][1];
			//beginning of integration block
			for(kk=kmin;kk<=kmax;kk++){
				xtotin[kk]=0.;
				xLPo[kk]=0.;
				xCP[kk]=0.;
				xEVPA[kk]=0.;
              //RG:FLAG SHOULDN'T IT BE:
				// totin[kk]=0.;
				// LPo[kk]=0.;
				// CP[kk]=0.;
				// EVPA[kk]=0.;
			};
			for(fnum=fmin;fnum<=fmax;fnum+=sep){
				init(sp,fmin,fmax,sep);           //initialization + parallel evaluation of intensity
				#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) shared(ittot)
				#include "intensity.cpp"
				for(kk=kmin;kk<=kmax;kk++){

					xtotin[kk]+=totin[kk]/ind;
					xLPo[kk]+=log(LPo[kk])/ind; //RG:FLAG WHY log() ? See Shcherbakov phd thesis
					xCP[kk]+=CP[kk]/ind;
					xEVPA[kk]+=EVPA[kk]/ind;

                    printf(YELLOW"[m_space.cpp]: "CYAN"oo=%d: xtotin[%d]=%f,totin[%d]=%f,ind=%d\n"RESET,oo,kk,xtotin[kk],kk,totin[kk],ind);

				};
			};
			//end of integration block 
			for(il=4;il<=10;il++)resid[il-4][oo]=(xtotin[il]-tofit[il][1])/dFnu[il]; //compute residuals with 1/dFnu weights
		};
		ans=(clock() - start ) / (doub)CLOCKS_PER_SEC;                               //timing
		printf(YELLOW"[m_space.cpp]: "RESET"Time = %.2f s; finished %d iterations; th=%.3f; heat=%.3f\n", ans,niter,th,heat);

		//evaluation of Jacobian, RHS and LHS of \chi^2 minimization
		int ij,ii;
		for(il=0;il<15;il++){Jac[il][0]=(resid[il][1]-resid[il][0])/drho;};
		for(ij=0;ij<1;ij++)for(ii=0;ii<1;ii++){matr[ij][ii]=0.;for(il=0;il<15;il++)matr[ij][ii]+=Jac[il][ij]*Jac[il][ii];};
		for(ii=0;ii<1;ii++){bb[ii]=0.;for(il=0;il<15;il++)bb[ii]+=Jac[il][ii]*resid[il][0];};
		ddr=-bb[0]/matr[0][0];                      //next relative offset of density
		printf(YELLOW"[m_space.cpp]: "RESET"Iterated ddr=%.4f\n",ddr);
		heat=inp[0][0];rhonor=inp[0][1]*(1+0.7*ddr);//due to high chance of overshooting and going to unphysical densities, we update density with 0.7*ddr - this slows the convergence, but makes it reliable!
	};

	doub xisq=0.,                      //\chi^2
		dof=7.;                        //degrees of freedom

    //RG:FLAG OK also when kmin!=0, but MUST HAVE kmin<5 !!!
	for(il=4;il<=10;il++)              //compute \chi^2
		xisq+=pow((doub)(xtotin[il]-tofit[il][1])/dFnu[il],(doub)2.);
	xisq/=dof;                         //reduced \chi^2
	stringstream sss;                  //filename for writing the "best-fitting" model parameters & spectra to disk
	sss<<(int)100*a<<"th"<<(int)floor(100*th+1e-6)<<"in"<<ind;
	string outstr=dir+"xesa"+sss.str(); 
	outstr+=".dat";
	
	FILE * xFile; 
	xFile=fopen(outstr.c_str(),"a");   //writing the "best-fitting" model parameters & spectra to disk
	if (heat==heat_min) fprintf(xFile,"# a \t th \t xisq \t 0.[?] \t heat \t rhonor \t rate*year/Msun \n");
	fprintf(xFile,"%.2f \t %.3f \t %.2f \t %.1f \t %.5f  %.5f \t %.3e\n",a,th,xisq,0.,heat,rhonor,rate*year/Msun);
	fclose(xFile);
	doub span=1.12;                    //ratio of consecutive density normalizations
	iswrite=true;                      //write output for a range of densities to disk

	rhonor/=span*span*span*span;       //vary density from rnohor/span^3 to rnohor*span^3
	for(w=0;w<7;w++){
		rhonor*=span;
		//beginning of spectrum calculation block
		for(kk=kmin;kk<=kmax;kk++){
			xtotin[kk]=0.;
			xLPo[kk]=0.;
			xCP[kk]=0.;
			xEVPA[kk]=0.;
		};
		for(fnum=fmin;fnum<=fmax;fnum+=sep){
			init(sp,fmin,fmax,sep);
			#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) shared(ittot)
			#include "intensity.cpp"
			for(kk=kmin;kk<=kmax;kk++){
				xtotin[kk]+=totin[kk]/ind;
				xLPo[kk]+=LPo[kk]/ind;
				xCP[kk]+=CP[kk]/ind;
				xEVPA[kk]+=EVPA[kk]/ind;
			};
		};
		//end of block
		//writing spectra into files and on screen
		string stra = sss.str();
		FILE * pFile;
		pFile = fopen ((dir+"ava"+stra+".dat").c_str(),"a");
        fprintf(pFile,"# sftab[kk][0]\t xtotin[kk]\txLPo[kk] xCP[kk] xEVPA[kk] heat \t rhonor Bpo\n");
		for(kk=kmin;kk<=kmax;kk++){
			printf(YELLOW"[m_space.cpp]: "RESET"avg at f=%.1f; I=%.3fJy LP=%.2f%% CP=%.3f%% EVPA=%.2fdeg\n",sftab[kk][0], xtotin[kk],xLPo[kk], xCP[kk],xEVPA[kk]);
			fprintf(pFile,"%.1f\t\t %.3f\t\t %.3f\t %.3f\t %.2f\t %.5f\t %.4f\n",sftab[kk][0], xtotin[kk],xLPo[kk], xCP[kk],xEVPA[kk],heat,rhonor,Bpo);
		};
		fclose(pFile);
	};
	rhonor/=span*span*span;//return rhonor to the best-fitting value
};
}
