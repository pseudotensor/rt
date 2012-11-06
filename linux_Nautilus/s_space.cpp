{//search for min xisq in "rhonor" for given "heat"
//jobs 1-5; sp=0..
	doub xaccur,xaccurr,xfact,xss,xsnxy,xstep,xsstep,xIint,xIang,xans;
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

cas=atoi(descr);printf("thx=%d\n",cas);//"LSB_JOBINDEX"
//if(argc==5){sp=atoi(argv[1])-1;co=atoi(argv[2]);if(atoi(argv[3])==0)boosted=false; else boosted=true;} else return 0;
printf("fnum+=%d; sp=%d\n",co,sp);//for various spins co should be set co=0;

	th=PI/2.*(2*cas-1)/thn;tth=th;xth=th;rhonor=3e5;//rhonor=3e5;// - for heat=0.7;
	for(heat=0.75;heat>0.15;heat*=0.96){
	//for(heat=0.35;heat>=0.35;heat*=0.96){
		ddr=1.;ddh=0.02;niter=0;iswrite=false;wfit[0]=false;wfit[1]=false;wfit[2]=false;
		while(fabs(ddr)>0.01){niter++;
		inp[0][0]=heat;inp[0][1]=rhonor;
		inp[1][0]=heat;inp[1][1]=rhonor*(1+drho);start=std::clock();

		for(oo=0;oo<2;oo++){
		heat=inp[oo][0];rhonor=inp[oo][1];
		for(fnum=-1-co;fnum<=-1/*8450*/-co;fnum+=sep){init(sp);
		#pragma omp parallel for schedule(dynamic,1) shared(ittot)
		#include "intensity.cpp"
		;};}
		ans=( std::clock() - start ) / (doub)CLOCKS_PER_SEC;printf ("Time = %.2f s; finished %d iterations; th=%.3f; heat=%.3f\n", ans,niter,th,heat);
	
		int ij,ii;//if(xxisq<lastxisq){lastxisq=xxisq;lastrhonor=inp[0][1];};
		for(il=0;il<15;il++){Jac[il][0]=(resid[il][1]-resid[il][0])/drho;};
		for(ij=0;ij<1;ij++)for(ii=0;ii<1;ii++){matr[ij][ii]=0.;for(il=0;il<15;il++)matr[ij][ii]+=Jac[il][ij]*Jac[il][ii];};
		for(ii=0;ii<1;ii++){bb[ii]=0.;for(il=0;il<15;il++)bb[ii]+=Jac[il][ii]*resid[il][0];};
		ddr=-bb[0]/matr[0][0];
		printf("Iterated ddr=%.4f\n",ddr);
		heat=inp[0][0];rhonor=inp[0][1]*(1+0.7*ddr);//remember 0.7!!!
		;};

		fnum=-1-co;oo=0;init(sp);//rhonor=lastrhonor;
		stringstream sss;sss<<(int)100*a<<"th"<<(int)floor(100*tth+1e-6);sss<<"fn"<<fnum;if(boosted)sss<<"boo"; else sss<<"nob";//to concide with last digit of fnum
		outstr=dir+"resa"+sss.str(); outstr+=".dat";
		FILE * xFile; xFile = fopen (outstr.c_str(),"a");
		fprintf(xFile,"%.2f %.3f %.5f %.5f %.5f %.5f %.3e\n",a,tth,xxisq,xang,heat,rhonor,rate*year/Msun);
		fclose(xFile);doub span=1.12;
	iswrite=true;wfit[0]=false;wfit[1]=false;wfit[2]=true;
	rhonor/=span*span*span*span;
	for(w=0;w<7;w++){
	rhonor*=span;
	for(fnum=-1-co;fnum<=-1-co;fnum+=sep){init(sp);
	#pragma omp parallel for schedule(dynamic,1) shared(ittot)
	#include "intensity.cpp"
	;};};rhonor/=span*span*span;//scanning the parameter space
/*
	sep=1;iswrite=false;averaged=false;kmin=co;kmax=co;
	for(fnum=6950;fnum<=9950;fnum+=sep){init(sp);
	#pragma omp parallel for schedule(dynamic,1) shared(ittot)
	#include "imaging.cpp"
	//#include "intensity.cpp"
*/
/*for(kk=kmin;kk<=kmax;kk++){
stringstream sstr;sstr <<(int)100*a<<"th"<<(int)floor(100*tth+1e-6)<<"nu"<<(int)floor(sftab[kk][0]+1e-6)<<"fn0";string stra = sstr.str();
FILE * pFile; pFile = fopen ((dir+"poliresa"+stra+".dat").c_str(),"a");
fprintf(pFile,"%d %.2f %.5f %.3f %.3f %.3f %.3f %.5f %.4f %.4f %.2f %.4e %.4f\n",
		fnum,sftab[kk][0],totin[kk],LPo[kk],ang[kk],CP[kk],err[kk],heat,rhonor,xxisq,TpTe,rate*year/Msun,tth);
fclose(pFile);}
}
*/

	};
}