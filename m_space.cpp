{//surf the entire parameter space searching for the best fit to polarized spectrum

  doub 
    xaccur=3e-3,  //1. absolute accuracy of geodesics computation
    xaccurr=1e-2, //2. absolute accuracy of radiative transfer integration
    xfact=1.0,    //3. relative size of the integration region
    xss=1e-2,     //4. fractional distance from BH horizon to the sphere, where geodesic integration stops
    xsnxy=101,    //5. number of points N along one side in the picture plane for N x N intensity grid
    xstep=0.1,    //6. step size in geodesic computation
    xsstep=-0.09, //7. step size in radiative transfer computation
    xIint=2e-9,   //8. initial intensity along each ray for radiative transfer
    xIang=1.,     //9. initial polarized phases along each ray for radiative transfer
    step=xstep,   //local variables, which control radiative transfer
    sstep=xsstep, 
    Iint=xIint,
    Iang=xIang,
    resid[sflen][4], //residuals for \chi^2 minimization algorithm
    Jac[sflen][3],   //correspondent Jacobian
    bb[3],           //RHS vector for the next minimization iteration
    matr[3][3];      //LHS matrix for the next minimization iteration

  int 
    snxy=xsnxy,    //global variable correspondent to xsnxy
 	ind=co,        //number of snapshots over which to compute mean intensity = 2-nd command line argument
 	oo,            // model label used together with inp array ~> [ASTRORAY_main.cpp]
    fmin=2150,     //minimum ID of fluid simulation snapshot
 	fmax=2150;     //maximum ID of fluid simulation snapshot
 
 // sep =1;

accur=xaccur;      //assigning values to global variables, which control radiative transfer
accurr=xaccurr;
fact=xfact;
ss=xss;

cas=atoi(descr);      //batch job ID = environment variable

//RG: DEFAULT: ASTRORAYv1
//th=PI/2.*(2*cas-1)/thn;//define theta angle out of a discrete set

start=clock();

rhonor=3e5;           //initial density for iterations
drho=0.02;            //density relative offset to compute Jacobian

switch (cas){
#include "lightup_jet.cpp"
}

sep=(fmax-fmin)/(ind-1);//ID difference between consecutive considered fluid simulation snapshots

printf(YELLOW"[m_space.cpp]: "RESET"in=%d shots; sp=%d; angle th=%f\n",co,sp,th);



/**************************************************************
 * STEP 1: Loop through heat and find best rhonor ~> best-fit *
 **************************************************************/


doub heat_min=0.3, heat_max=0.75;
// for(heat=heat_max; heat>heat_min; heat*=0.95){     // choose a set of heat parameter
for(heat=heat_min; heat<=heat_max; heat+=0.1){     // choose a set of heat parameter

	doub ddr=1.;                                   // relative adjustment of density for next iteration
	niter=0;                                       // number of iterations
	iswrite=false;
	while(fabs(ddr)>0.01){                         // do iterations till convergence is reached
		niter++;

        // initialize 2 states with same heat and slightly different rhonor
		inp[0][0]=heat;                          
		inp[0][1]=rhonor;
		inp[1][0]=heat;
		inp[1][1]=rhonor*(1+drho);

		// compute spectra for both states
		for(oo=0;oo<2;oo++){                      
          heat=inp[oo][0]; // see [ASTRORAY_main.cpp] inp[][]: 1st dim stores different models; 2nd dim stores free parameters: heat,rhonor,th
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


            // case 42351: fmin=2140;fmax=2160;sep=1;kmin=4;kmax=10;sp=0;rhonor=4.5e7; heat=0.46; th=PI/2-0.4; dphi=4.*PI/3.;thlimit=0.05;fdiff=40;isBcut=false;isBred=true;magn_cap=10.;Te_jet_par=35.;break;
            // fmin=2140;fmax=2160;sep=1;kmin=4;kmax=10;sp=0;rhonor=4.5e7; heat=0.46; th=PI/2-0.4; dphi=4.*PI/3.;thlimit=0.05;fdiff=40;isBcut=false;isBred=true;magn_cap=10.;Te_jet_par=35.;


			for(fnum=fmin;fnum<=fmax;fnum+=sep){

                // ID; QUERY GRMHD DATA
				init(sp,fmin,fmax,sep);

                // parallel evaluation of intensity
				#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) shared(ittot)
				#include "intensity.cpp"

                // Take averages over different snapshots
				for(kk=kmin;kk<=kmax;kk++){

					xtotin[kk]+=totin[kk]/ind;
					xLPo[kk]+=log(LPo[kk])/ind; //RG:FLAG WHY log()? Force it to be positive? See Shcherbakov phd thesis...
					xCP[kk]+=CP[kk]/ind;
					xEVPA[kk]+=EVPA[kk]/ind;

                    printf(YELLOW"[m_space.cpp]: "CYAN"oo=%d: xtotin[%d]=%f,totin[%d]=%f,ind=%d\n"RESET,oo,kk,xtotin[kk],kk,totin[kk],ind);

				}; // for(kk=kmin;kk<=kmax;kk++){
			}; // for(fnum=fmin;fnum<=fmax;fnum+=sep){



            // RESIDUALS

			for(il=4;il<=10;il++)resid[il-4][oo]=(xtotin[il]-tofit[il][1])/dFnu[il]; //compute residuals with 1/dFnu weights
		};

		ans=(clock() - start ) / (doub)CLOCKS_PER_SEC;                               //timing
		printf(YELLOW"[m_space.cpp]: "RESET"Time = %.2f s; finished %d iterations; th=%.3f; heat=%.3f\n", ans,niter,th,heat); //report

		//evaluation of Jacobian, RHS and LHS of \chi^2 minimization
		int ij,ii;
        //RG:FLAG WHY 15?
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
	sss<<(int)100*a<<"th"<<(int)floor(100*th+1e-6)<<"in"<<ind<<"case"<<cas;
	string outstr=dir+"xesa"+sss.str(); 
	outstr+=".dat";
	
	FILE * xFile; 
	xFile=fopen(outstr.c_str(),"a");   //writing the "best-fitting" model parameters & spectra to disk
	if (heat==heat_min) fprintf(xFile,"# a \t th \t xisq \t 0.[?] \t heat \t rhonor \t rate*year/Msun \n");
	fprintf(xFile,"%.2f \t %.3f \t %.2f \t %.1f \t %.5f %.5f  %.3e\n",a,th,xisq,0.,heat,rhonor,rate*year/Msun);
	fclose(xFile);



    /******************************************************************
     * STEP 2: Vary rhonor around best-fit obtained in STEP 1 (above) *
     ******************************************************************/


	// doub rho_scan_factor=1.12;                    //ratio of consecutive density normalizations
	doub rho_scan_factor=1.1;                    // steps/increments in density normalizations to scan for
	iswrite=true;                      //write output for a range of densities to disk

    //vary density between [rhonor/rho_scan_factor^3 , rhonor*rho_scan_factor^3]
	// rhonor/=rho_scan_factor*rho_scan_factor*rho_scan_factor*rho_scan_factor;
    doub rho_scan_min=pow(rhonor,-3*rho_scan_factor);
    doub rho_scan_max=pow(rhonor,+3*rho_scan_factor);

    //RG:FLAG WHAT IS "w"? hardcoded nr of iterations in rhonor scan?! Does not appear anywhere after...
	// for(w=0;w<7;w++){ // in [ASTRORAY_main.cpp]: "int w; //thread number (for testing)"
	// 	rhonor*=rho_scan_factor;

    for(doub rho_scan=rho_scan_min; rho_scan>=rho_scan_max; rho_scan*=rho_scan_factor){ // Why not this way??
        rhonor=rho_scan;
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
        fprintf(pFile,"# nu\t <I>\t<LP> <CP> <EVPA> heat \t rhonor Bpo\n");
		for(kk=kmin;kk<=kmax;kk++){
			printf(YELLOW"[m_space.cpp]: "RESET"avg at f=%.1f; I=%.3fJy LP=%.2f%% CP=%.3f%% EVPA=%.2fdeg\n",sftab[kk][0], xtotin[kk],xLPo[kk], xCP[kk],xEVPA[kk]);
			fprintf(pFile,"%.1f\t\t %.3f\t\t %.3f\t %.3f\t %.2f\t %.5f\t %.4f\n",sftab[kk][0], xtotin[kk],xLPo[kk], xCP[kk],xEVPA[kk],heat,rhonor,Bpo);
		};
		fclose(pFile);
	};

    //RG:FLAG WHY DIVIDE *THREE TIMES* BY "rho_scan_factor" ?
	rhonor/=rho_scan_factor*rho_scan_factor*rho_scan_factor;//return rhonor to the best-fitting value
};
}
