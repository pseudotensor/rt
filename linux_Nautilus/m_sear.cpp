{doub xaccur,xaccurr,xfact,xss,xsnxy,xstep,xsstep,xIint,xIang,xans;int testN, ind;
xaccur=3e-3; //+1. accuracy of geodesics computation
xaccurr=1e-2;//+2. accuracy of radiative transfer integration
xfact=1.0;   //3. size of the integration region
xss=1e-2;    //+4. fractional distance from BH horizon to start integration
xsnxy=111;    //+5. grid N in the picture plane
xstep=0.1;   //+6. step control of geodesic computation
xsstep=-0.09;//+7. step control of radiative transfer computation
xIint=2e-9;  //+8. initial intensity along each ray
xIang=1.;    //+9. initial polarized phases along each ray

//xsnxy=55;//for testing

accur=xaccur;accurr=xaccurr;fact=xfact;ss=xss;int snxy=xsnxy;doub step=xstep,sstep=xsstep,Iint=xIint,Iang=xIang,devsq=0.;


int fmin,fmax,sep;
//ind=21;fdiff=60;
fmin=5534;fmax=7034;sp=0;rhonor=144694.13267; heat=0.17339;th=2.4059;thlimit=0.0;fdiff=0;isBcut=true;//Jon simulations with B^2/\rho cut-off
fmin=5534;fmax=7034;sp=0;rhonor=192642.42683 ; heat=0.16737;th=2.61;thlimit=0.1;fdiff=0;isBcut=true;//Jon simulations, B^2/\rho cut-off, thlimit=0.1

fmin=5534;fmax=7034;sp=0;rhonor=142642.42683;heat=0.16737;th=2.6367;thlimit=0.1;fdiff=0;isBcut=false;isBred=false;//Jon simulations, no B^2/\rho cut-off, thlimit=0.1
fmin=5534;fmax=7034;sp=0;rhonor=288331.85590;heat=0.14320;th=2.4529;thlimit=0.1;fdiff=0;isBcut=false;isBred=true;//Jon simulations, B^2/\rho reduction, thlimit=0.1, no B^2\rho cutoff

fmin=2000;fmax=7000;sp=0;rhonor=261385.84479;heat=0.15572;th=2.3768;thlimit=0.1;fdiff=(int)Bpo; isBcut=false;isBred=true;//Jon simulations, B^2/\rho reduction, thlimit=0.1, no B^2\rho cutoff
fmin=2000;fmax=7000;sp=0;rhonor=144897.83871;heat=0.17172;th=2.3897;thlimit=0.1;fdiff=(int)Bpo; isBcut=false;isBred=false;//Jon simulations, thlimit=0.1, no B^2\rho cutoff or B^2/\rho reduction
fmin=2000;fmax=7000;sp=0;rhonor=144897.83871;heat=0.17172;th=2.3897;thlimit=0.1;fdiff=(int)Bpo; isBcut=false;isBred=false;//Olek simulations thlimit=0.1;

//if(co==21){fmin=6850;fmax=9850;}
//else{
//temporary block for \chi^2 over various time windows
//switch (co){
//case 1: fmin=6850;fmax=7350;break;											
//case 2: fmin=7350;fmax=7850;break;								
//case 3: fmin=7850;fmax=8350;break;								
//case 4:	fmin=8350;fmax=8850;break;						
//case 5:	fmin=8850;fmax=9350;break;						
//case 6:	fmin=9350;fmax=9850;break;
//}
//}
//end of temporary block

//sep=4000;ind=1;//for testing sep>fmax-fmin

ind=co;testN=atoi(descr);
nP=12;if(isLP87)nPeff=nP; else nPeff=nP-1;
sep=(fmax-fmin)/(ind-1);
printf("testN=%d;  ind=%d\n",testN,ind);

dheat=0.;drho=0.;dtheta=0.;
switch (testN){
case 1: break;											
case 2: dheat=0.15;dheat=0.9;break;								
case 3: dheat=-0.15;dheat=-0.5;break;								
case 4:	drho=0.25;drho=0.9;break;						
case 5:	drho=-0.25;drho=-0.6;break;
case 6:	dtheta=0.1;dtheta=0.5;break;
case 7:	dtheta=-0.1;dtheta=-0.5;break;
}
/*switch (testN){
case 1: th=2.4840;heat=0.16992;rhonor=147780.66084;break;
case 2: th=2.4812;heat=0.21682;rhonor=104110.02380;break;
case 3: th=2.4198;heat=0.13443;rhonor=171629.22288;break;
case 4: th=2.4978;heat=0.16175;rhonor=163531.22261;break;
case 5: th=2.4039;heat=0.17473;rhonor=140163.50164;break;
case 6: th=2.6143;heat=0.16798;rhonor=149993.46045;break;
case 7: th=2.1688;heat=0.16663;rhonor=156507.65524;break;
}*/



heat*=(1+dheat);rhonor*=(1+drho);th+=dtheta;
	//sp=5;rhonor=361204.; heat=0.445612; th=2.125;
printf("in=%d shots; sp=%d; th=%.4f; heat=%.4f; rhonor=%.1f; fdiff=%d\n",co,sp,th,heat,rhonor,fdiff);//


		ddr=1.;ddh=0.02;dth=0.01;niter=0;iswrite=false;lastxisq=8.;tth=th;xth=th;
		dheat=0.02;drho=0.02;dtheta=0.01;// based on convergence tests
while((fabs(ddh)>0.003)||(fabs(ddr)>0.01)||(fabs(dth)>0.005)){niter++;
		inp[0][0]=heat;inp[0][1]=rhonor;inp[0][2]=th;
		inp[1][0]=heat*(1+dheat);inp[1][1]=rhonor;inp[1][2]=th;
		inp[2][0]=heat;inp[2][1]=rhonor*(1+drho);inp[2][2]=th;
		inp[3][0]=heat;inp[3][1]=rhonor;inp[3][2]=th*(1+dtheta);
		start=std::clock();
	 for(oo=0;oo<4;oo++)for(kk=kmin;kk<=kmax;kk++){ztotin[kk][oo]=0.;zLPo[kk][oo]=0.;zCP[kk][oo]=0.;zEVPA[kk][oo]=0.;};
//beginning of integration block
		for(fnum=fmin;fnum<=fmax;fnum+=sep){
			for(oo=0;oo<4;oo++){//2D heat and rhonor integration converges quickly, whereas 3D heat, rhonor and th integration has problems (since LP changes unpredictably with changes in th)
		heat=inp[oo][0];rhonor=inp[oo][1];th=inp[oo][2];
		init(sp,fmin,fmax,sep);
		#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) shared(ittot)
		#include "intensity.cpp"
		for(kk=kmin;kk<=kmax;kk++){ztotin[kk][oo]+=totin[kk]/ind;zLPo[kk][oo]+=LPo[kk]/ind;zCP[kk][oo]+=CP[kk]/ind;zEVPA[kk][oo]+=ang[kk]/ind;};
		;};}
//end of integration block 
	for(oo=0;oo<4;oo++){
		for(il=4;il<=10;il++)resid[il-4][oo]=(ztotin[il][oo]-tofit[il][1])/dFnu[il];
		resid[7][oo]=(zLPo[4][oo]-tofit[4][2])/dLP[0];
		resid[8][oo]=(zLPo[7][oo]-tofit[7][2])/dLP[1];
		resid[9][oo]=(zLPo[8][oo]-tofit[8][2])/dLP[2];
		resid[10][oo]=(zCP[7][oo]-tofit[7][4])/dCP;
		resid[11][oo]=(zCP[8][oo]-tofit[8][4])/dCP;}
{xisq=0.;for(il=0;il<nP;il++)xisq+=resid[il][0]*resid[il][0]/(nPeff-3);
			printf("heat=%.4f, rhonor=%.1f, th=%.4f, xisq=%.3f\n",inp[0][0], inp[0][1], inp[0][2], xisq);
			stringstream sss;sss<<(int)100*a<<"in"<<co<<"N"<<testN<<"fdiff"<<fdiff;//"th"<<(int)floor(100*tth+1e-6)<<
			outstr=dir+"xisqa"+sss.str(); outstr+=".dat";
			FILE * xFile; xFile = fopen (outstr.c_str(),"a");
			fprintf(xFile,"%.2f %.4f %.4f %.5f %.5f %.3e\n",a,xisq,inp[0][2],inp[0][0],inp[0][1],rate*year/Msun);fclose(xFile);};

		ans=( std::clock() - start ) / (doub)CLOCKS_PER_SEC;printf ("Time = %.2f s; finished %d iterations; th=%.3f; heat=%.3f\n", ans,niter,th,heat);
	

		int ij,ii;//if(xxisq<lastxisq){lastxisq=xxisq;lastrhonor=inp[0][1];lastheat=inp[0][0];};
		for(il=0;il<nP;il++){Jac[il][0]=(resid[il][1]-resid[il][0])/dheat;Jac[il][1]=(resid[il][2]-resid[il][0])/drho;Jac[il][2]=(resid[il][3]-resid[il][0])/dtheta;};
		for(ij=0;ij<3;ij++)for(ii=0;ii<3;ii++){matr[ij][ii]=0.;for(il=0;il<nP;il++)matr[ij][ii]+=Jac[il][ij]*Jac[il][ii];};
		for(il=0;il<3;il++)matr[il][il]*=1.1;//Levenberg-Marquardt

		for(ii=0;ii<3;ii++){bb[ii]=0.;for(il=0;il<nP;il++)bb[ii]+=Jac[il][ii]*resid[il][0];};
det=-(matr[0][2]*matr[1][1]*matr[2][0]) + matr[0][1]*matr[1][2]*matr[2][0] + matr[0][2]*matr[1][0]*matr[2][1] - matr[0][0]*matr[1][2]*matr[2][1] - matr[0][1]*matr[1][0]*matr[2][2] + matr[0][0]*matr[1][1]*matr[2][2];
ddh=(bb[2]*matr[0][2]*matr[1][1] - bb[2]*matr[0][1]*matr[1][2] - bb[1]*matr[0][2]*matr[2][1] + bb[0]*matr[1][2]*matr[2][1] + (bb[1]*matr[0][1] - bb[0]*matr[1][1])*matr[2][2])/det;
ddr=(bb[2]*(-(matr[0][2]*matr[1][0]) + matr[0][0]*matr[1][2]) + bb[1]*(matr[0][2]*matr[2][0] - matr[0][0]*matr[2][2]) + bb[0]*(-(matr[1][2]*matr[2][0]) + matr[1][0]*matr[2][2]))/det;
dth=(bb[2]*matr[0][1]*matr[1][0] - bb[2]*matr[0][0]*matr[1][1] - bb[1]*matr[0][1]*matr[2][0] + bb[0]*matr[1][1]*matr[2][0] + (bb[1]*matr[0][0] - bb[0]*matr[1][0])*matr[2][1])/det;
		//ddh=-((-(bb[1]*matr[0][1]) + bb[0]*matr[1][1])/(-matr[0][1]*matr[0][1] + matr[0][0]*matr[1][1]));
		//ddr=-((bb[1]*matr[0][0] - bb[0]*matr[0][1])/(-matr[0][1]*matr[0][1] + matr[0][0]*matr[1][1]));
		printf("Iterated ddr=%.4f, ddh=%.4f, dth=%.4f\n",ddr,ddh,dth);
		while((fabs(ddh)>0.10) || (fabs(ddr)>0.2) || (fabs(dth)>0.05) /*|| (inp[0][2]*(1+det*dth)>PI) || (inp[0][2]*(1+det*dth)<0)*/){ddh/=2.;ddr/=2.;dth/=2.;}
	det=1.0; //remember 0.7!!!
		heat=inp[0][0]*(1+det*ddh);rhonor=inp[0][1]*(1+det*ddr);th=inp[0][2]*(1+det*dth);
		if (niter>20){break;}
;};

stringstream sss;sss<<(int)100*a;//<<"th"<<(int)floor(100*tth+1e-6)<<"in"<<ind

	//beginning of integration block
	for(kk=kmin;kk<=kmax;kk++){xtotin[kk]=0.;xLPo[kk]=0.;xCP[kk]=0.;xEVPA[kk]=0.;};
	for(fnum=fmin;fnum<=fmax;fnum+=sep){init(sp,fmin,fmax,sep);
	#pragma omp parallel for schedule(dynamic,1) num_threads(nthreads) shared(ittot)
	#include "intensity.cpp"
	for(kk=kmin;kk<=kmax;kk++){xtotin[kk]+=totin[kk]/ind;xLPo[kk]+=LPo[kk]/ind;xCP[kk]+=CP[kk]/ind;xEVPA[kk]+=ang[kk]/ind;};
	;}
	//end of integration block 	
		for(il=4;il<=10;il++)resid[il-4][0]=(xtotin[il]-tofit[il][1])/dFnu[il];
		resid[7][oo]=(xLPo[4]-tofit[4][2])/dLP[0];
		resid[8][oo]=(xLPo[7]-tofit[7][2])/dLP[1];
		resid[9][oo]=(xLPo[8]-tofit[8][2])/dLP[2];
		resid[10][oo]=(xCP[7]-tofit[7][4])/dCP;
		resid[11][oo]=(xCP[8]-tofit[8][4])/dCP;
		xisq=0.;for(il=0;il<nP;il++)xisq+=resid[il][0]*resid[il][0]/(nPeff-3);

string stra = sss.str();FILE * pFile; 
pFile = fopen ((dir+"bestfita"+stra+".dat").c_str(),"a");
	for(kk=kmin;kk<=kmax;kk++){
		printf("avg at f=%.1f; I=%.3fJy LP=%.2f%% CP=%.3f%% EVPA=%.2fdeg\n",sftab[kk][0], xtotin[kk],xLPo[kk], xCP[kk],xEVPA[kk]);
		fprintf(pFile,"%.1f %.3f %.3f %.3f %.2f %.5f %.4f %.4f %.4f\n",sftab[kk][0], xtotin[kk],xLPo[kk], xCP[kk],xEVPA[kk],th,heat,rhonor,Bpo);
		};fclose(pFile);

}

/* ---- found best \chi^2/dof
Jon's simulations thickdisk7:
no manipulations        ; 6434-7034 range in  2: 0.94 6.0003 2.1166 0.18630 93055.86347 8.355e-009 //1
no manipulations        ; 6434-7034 range in 21: 0.94 7.5319 1.8225 0.12172 210453.42461 1.890e-008//1
B^2/rho precise cutoff  ; 5534-7034 range in  2: 0.94 4.8945 2.1361 0.14852 250637.00991 2.250e-008//2
costh=0.05 polar cutoff ; 5534-7034 range in  2: 0.94 4.5160 2.2646 0.17137 140141.68941 1.258e-008//3
costh=0.1  polar cutoff ; 5534-7034 range in  2: 0.94 2.6097 2.4059 0.17339 144694.13267 1.299e-008//4
*/
/*
case 1: break;											//0.50 2.0997 2.0298 0.38745 912995.17514 2.678e-08
case 2: dheat=0.1;break;								//0.50 2.1435 2.0266 0.38579 919132.55441 2.696e-08
case 3: dheat=-0.2;break;								//0.50 2.2671 2.0345 0.37847 957820.38936 2.809e-08
case 4:	dheat=0.3;drho=0.2;break;						//0.50 2.1162 2.0274 0.38548 921203.11770 2.702e-08
case 5:	dheat=-0.3;drho=0.1;break;						//0.50 5.1026 2.0203 0.33018 1265744.41428 3.712e-08
case 6:	dheat=-0.3;drho=0.1;dtheta=0.05;break;//-++		//0.50 3.7341 2.1323 0.36096 1049217.99214 3.077e-08
case 7:	dheat=0.3;drho=-0.2;dtheta=-0.05;break;//+--	//0.50 2.6665 2.0452 0.41094 810315.66064 2.376e-08
case 8:	dheat=-0.4;drho=-0.1;dtheta=0.1;break;//--+		//0.50 5.0790 2.2129 0.33548 1229052.93333 3.604e-08
case 9: dheat=-0.3;drho=0.1;dtheta=-0.05;break;//-+-	//0.50 3.1288 1.9725 0.37329 976050.33515 2.862e-08
case 10: dheat=-0.3;drho=-0.1;dtheta=-0.05;break;//---	//0.50 2.9477 1.9595 0.37366 975031.84604 2.860e-08
case 11: dheat=0.3;drho=0.1;dtheta=0.08;break;//+++		//0.50 4.7750 2.1661 0.40325 855991.69689 2.510e-08//not final
case 12: dheat=0.3;drho=-0.15;dtheta=0.08;break;//+-+	//0.50 3.5087 2.0271 0.43631 733988.87914 2.153e-08
*/