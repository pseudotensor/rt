{//search for min xisq in "rhonor" for given "heat" //jobs 1-5; sp=0..
doub xaccur,xaccurr,xfact,xss,xsnxy,xstep,xsstep,xIint,xIang,xans;
xaccur=3e-3; //+1. accuracy of geodesics computation
xaccurr=1e-2;//+2. accuracy of radiative transfer integration
xfact=1.0;   //3. size of the integration region
xss=1e-2;    //+4. fractional distance from BH horizon to start integration
xsnxy=111;   //+5. grid N in the picture plane
xstep=0.1;   //+6. step control of geodesic computation
xsstep=-0.09;//+7. step control of radiative transfer computation
xIint=2e-9;  //+8. initial intensity along each ray
xIang=1.;    //+9. initial polarized phases along each ray
accur=xaccur;accurr=xaccurr;fact=xfact;ss=xss;int snxy=xsnxy;doub step=xstep,sstep=xsstep,Iint=xIint,Iang=xIang,devsq=0.;

cas=atoi(descr);printf("angle thx=%d\n",cas);//"LSB_JOBINDEX"
printf("in=%d shots; sp=%d\n",co,sp);//for various spins co should be set co=0;

int ind=co;
int fmin=6950,fmax=9950,sep=(fmax-fmin)/(ind-1);drho=0.02;
	th=PI/2.*(2*cas-1)/thn;tth=th;xth=th;rhonor=3e5;//rhonor=3e5;// - for heat=0.7;

	for(heat=0.75;heat>0.15;heat*=0.96){
	//for(heat=0.35;heat>=0.35;heat*=0.96){
		ddr=1.;ddh=0.02;niter=0;iswrite=false;
		while(fabs(ddr)>0.01){niter++;
		inp[0][0]=heat;inp[0][1]=rhonor;
		inp[1][0]=heat;inp[1][1]=rhonor*(1+drho);start=std::clock();

		for(oo=0;oo<2;oo++){
		heat=inp[oo][0];rhonor=inp[oo][1];
		//beginning of integration block
		for(kk=kmin;kk<=kmax;kk++){xtotin[kk]=0.;xLPo[kk]=0.;xCP[kk]=0.;xEVPA[kk]=0.;};
		for(fnum=fmin;fnum<=fmax;fnum+=sep){init(sp,fmin,fmax,sep);
		#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) shared(ittot)
		#include "intensity.cpp"
		for(kk=kmin;kk<=kmax;kk++){xtotin[kk]+=totin[kk]/ind;xLPo[kk]+=log(LPo[kk])/ind;xCP[kk]+=CP[kk]/ind;xEVPA[kk]+=ang[kk]/ind;};
		;}
		//end of integration block 
		for(il=4;il<=10;il++)resid[il-4][oo]=(xtotin[il]-tofit[il][1])/dFnu[il];
		;}
		ans=( std::clock() - start ) / (doub)CLOCKS_PER_SEC;printf ("Time = %.2f s; finished %d iterations; th=%.3f; heat=%.3f\n", ans,niter,th,heat);
	
		int ij,ii;
		for(il=0;il<15;il++){Jac[il][0]=(resid[il][1]-resid[il][0])/drho;};
		for(ij=0;ij<1;ij++)for(ii=0;ii<1;ii++){matr[ij][ii]=0.;for(il=0;il<15;il++)matr[ij][ii]+=Jac[il][ij]*Jac[il][ii];};
		for(ii=0;ii<1;ii++){bb[ii]=0.;for(il=0;il<15;il++)bb[ii]+=Jac[il][ii]*resid[il][0];};
		ddr=-bb[0]/matr[0][0];
		printf("Iterated ddr=%.4f\n",ddr);
		heat=inp[0][0];rhonor=inp[0][1]*(1+0.7*ddr);//remember 0.7!!!
		;};

		fnum=fmin;oo=0;init(sp,fmin,fmax,sep);//rhonor=lastrhonor;
xisq=0;for(il=4;il<=10;il++)xisq+=pow((doub)(xtotin[il]-tofit[il][1])/dFnu[il],(doub)2.);dof=7;xisq/=dof;
stringstream sss;sss<<(int)100*a<<"th"<<(int)floor(100*tth+1e-6)<<"in"<<ind;
outstr=dir+"xesa"+sss.str(); outstr+=".dat";
FILE * xFile; xFile=fopen(outstr.c_str(),"a");
		fprintf(xFile,"%.2f %.3f %.5f %.5f %.5f %.5f %.3e\n",a,tth,xisq,xang,heat,rhonor,rate*year/Msun);fclose(xFile);
		doub span=1.12;
	iswrite=true;
	rhonor/=span*span*span*span;
	for(w=0;w<7;w++){
	rhonor*=span;
	//beginning of integration block
	for(kk=kmin;kk<=kmax;kk++){xtotin[kk]=0.;xLPo[kk]=0.;xCP[kk]=0.;xEVPA[kk]=0.;};
	for(fnum=fmin;fnum<=fmax;fnum+=sep){init(sp,fmin,fmax,sep);
	#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) shared(ittot)
	#include "intensity.cpp"
	for(kk=kmin;kk<=kmax;kk++){xtotin[kk]+=totin[kk]/ind;xLPo[kk]+=LPo[kk]/ind;xCP[kk]+=CP[kk]/ind;xEVPA[kk]+=ang[kk]/ind;};
	;}
	//end of integration block 	
string stra = sss.str();FILE * pFile; 
pFile = fopen ((dir+"ava"+stra+".dat").c_str(),"a");
	for(kk=kmin;kk<=kmax;kk++){
		printf("avg at f=%.1f; I=%.3fJy LP=%.2f%% CP=%.3f%% EVPA=%.2fdeg\n",sftab[kk][0], xtotin[kk],xLPo[kk], xCP[kk],xEVPA[kk]);
		fprintf(pFile,"%.1f %.3f %.3f %.3f %.2f %.5f %.4f\n",sftab[kk][0], xtotin[kk],xLPo[kk], xCP[kk],xEVPA[kk],heat,rhonor,Bpo);
		};fclose(pFile);
	};rhonor/=span*span*span;//scanning the parameter space
	};
}