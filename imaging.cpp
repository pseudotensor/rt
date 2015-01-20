for(ix=0;ix<(nxy+1)*(nxy+1);ix++){                             //cycle over all geodesics (more efficient with OpenMP than cycling along each direction in a double for loop)

  //RG: 
  //printf("ix=%d",ix);

	    int kk,
		iiy=ix % (nxy+1),                                      //find geodesic coordinates along x and y directions
		iix=(int)ix/(snxy+1);
	for(kk=kmin;kk<=kmax;kk+=kstep){                           //cycle over chosen frequencies
		doub t,
			 maxy=fact*sftab[kk][1],                           //get a size in picture plane for each frequency
			 xg=-maxy+2.*maxy/nxy*doub(iix),                   //offset of current geodesic in x direction
			 yg=-maxy+2.*maxy/nxy*doub(iiy),                   //offset of current geodesic in y direction
			 b=sqrt(xg*xg+yg*yg),                              //impact parameter
			 beta=atan2(yg,xg);                                //polar angle in picture plane
		//new geodesic integration for each ray
		#include "geoint.cpp"
		ppy[currth].nu=1e9*sftab[kk][0];                       //frequency
		//new radiative transfer is solved for each ray
		#include "solvetrans.cpp"
		
		(*ausin)[iix][iiy][kk][0]=II[0];                       //total intensity I
		(*ausin)[iix][iiy][kk][1]=II[1]*cos(II[2])*sin(II[3]); //linearly polarized intensity Q
		(*ausin)[iix][iiy][kk][2]=II[1]*sin(II[2])*sin(II[3]); //linearly polarized intensity U
		(*ausin)[iix][iiy][kk][3]=II[1]*cos(II[3]);            //circularly polarized intensity V
		(*ausin)[iix][iiy][kk][4]=II[4];                       //any other quantity integrated over a geodesic, i.e., approximate total intensity
	};
};

typedef double (*arra)[nxy+1][nxy+1][5];                       //array of intensities
   arra intab = (arra) new double[nxy+1][nxy+1][5];
 typedef double (*para)[20];                                   //array of parameters
  para params = (para) new double[20];
doub in[sflen][5];                                             //define array for computing intensity

for(kk=kmin;kk<=kmax;kk+=kstep){
	for(ix=0;ix<=nxy;ix++)
		for(iy=0;iy<=nxy;iy++)
			for(il=0;il<=4;il++)
				(*intab)[ix][iy][il]=(*ausin)[ix][iy][kk][il]; //separate arrays for different frequencies
	for(il=0;il<5;il++)                                        //initialize array of total intensities with zeros
		in[kk][il]=0.;
	doub hei=2./nxy,                                           //distance between neighbooring points in picture plane
		 maxy=fact*sftab[kk][1];                               //size of the integration regions in picture plane 

	for(ix=0;ix<=nxy-3;ix+=2)
		for(iy=0;iy<=nxy-3;iy+=2)
			for(il=0;il<=4;il++)                               //for each kind of intensity
				//2nd order accuracy of the integrator - works best
				in[kk][il]+=((*intab)[ix][iy][il] + (*intab)[ix][2 + iy][il] + 4*(*intab)[1 +ix][iy][il] + 4*(*intab)[1 + ix][2 + iy][il] + (*intab)[2 +ix][iy][il] +
				4*((*intab)[ix][1 + iy][il] + 4*(*intab)[1 + ix][1 +iy][il] + (*intab)[2 + ix][1 + iy][il]) + (*intab)[2 + ix][2 +iy][il])*hei*hei/9.;
	for(ix=0;ix<=nxy-1;ix++)                                   //integration over x boundary - 1-st order
		for(il=0;il<5;il++)
			in[kk][il]+=((*intab)[ix][nxy-1][il] + (*intab)[ix][nxy][il] + (*intab)[1 + ix][nxy-1][il] + (*intab)[1 + ix][nxy][il])*hei*hei/4.;
	for(iy=0;iy<=nxy-2;iy++)                                   //integration over y boundary - 1-st order
		for(il=0;il<5;il++)
			in[kk][il]+=((*intab)[nxy-1][iy][il] + (*intab)[nxy-1][1 + iy][il] + (*intab)[nxy][iy][il] + (*intab)[nxy][1 + iy][il])*hei*hei/4.;

	for(il=0;il<5;il++)
		in[kk][il]*=maxy*maxy;                                //normalization over integration region size
	totin[kk]=66.4648*in[kk][0];                              //normalization for angular size of Sgr A* region - recompute for your object!
	LPo[kk]=100.*sqrt(in[kk][1]*in[kk][1]+in[kk][2]*in[kk][2])/in[kk][0];//compute total LP fraction
	CP[kk]=100.*in[kk][3]/in[kk][0];                                     //compute CP fraction
	EVPA[kk]=fmod(180/PI*atan2(in[kk][2],in[kk][1])/2.+180.,180);        //compute EVPA
	err[kk]=66.4648*in[kk][4];                                           //normalize 5-th intensity as total intensity
	printf("%d; f=%.1f; I=%.3fJy LP=%.2f%% EVPA=%.1fdeg CP=%.3f%% non-pol I=%.2fJy \n",fnum,sftab[kk][0], totin[kk],LPo[kk],EVPA[kk], CP[kk],err[kk]);

	//setting parameters to be written in file
	(*params)[0]=a;
	(*params)[1]=th;
	(*params)[2]=double(nxy);
	(*params)[3]=double(sftab[kk][0]);
	(*params)[4]=double(sftab[kk][1]);
	(*params)[5]=double(heat);
	(*params)[6]=double(rhonor);
	(*params)[7]=double(totin[kk]);
	(*params)[8]=double(LPo[kk]);
	(*params)[9]=double(EVPA[kk]);
	(*params)[10]=double(CP[kk]);
	(*params)[11]=double(err[kk]);
	(*params)[12]=double(TpTe);
	(*params)[13]=rate*year/Msun;

	stringstream sstr;                                                             //prepare to write images into "shotimag" files
	sstr<<(int)100*a<<"th"<<(int)floor(100*th+1e-6)<<"f"<<floor(sftab[kk][0])<<"fn"<<fnum<<"_"<<nxy;
	string stra = sstr.str();

	ofstream faire ((dir+"shotimag"+stra+".dat").c_str(), ios::out|ios::binary);   //write into binary output file
	faire.write(reinterpret_cast<char *>(params), 20*sizeof(double));              //write params - fixed size variable
	faire.write(reinterpret_cast<char *>(intab), 5*(nxy+1)*(nxy+1)*sizeof(doub));  //write image
	faire.close();
}

if(iswrite){
	stringstream sstr;                                                             //prepare to write full spectra into "poliresa" files
	sstr <<(int)100*a<<"th"<<(int)floor(100*th+1e-6)<<"fn"<<fnum<<"hi";//"hi" is a for high resolution
	string stra = sstr.str();
	FILE * pFile; 
	pFile=fopen ((dir+"poliresa"+stra+fif+".dat").c_str(),"a");                    //open polires*
	for(kk=kmin;kk<=kmax;kk+=kstep)                                                //write full spectra into "poliresa" files
		fprintf(pFile,"%d %.2f %.5f %.4f %.4f %.4f %.4f %.5f %.4f %.4f %.2f %.4e %.4f\n", fnum,sftab[kk][0],totin[kk],LPo[kk],EVPA[kk],CP[kk],err[kk],heat,rhonor,0.,TpTe,rate*year/Msun,th);
	fclose(pFile);
};
delete [] intab;                                                                   //delete dynamic variable
delete [] params;                                                                  //delete parameters
