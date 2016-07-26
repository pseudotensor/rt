for(ix=0;ix<(nxy+1)*(nxy+1);ix++){                             //cycle over all geodesics (more efficient with OpenMP than cycling along each direction in a double for loop)

	    int kk,
		iiy=ix % (nxy+1),                                      //find geodesic coordinates along x and y directions
		iix=(int)ix/(nxy+1);
	for(kk=kmin;kk<=kmax;kk+=kstep){                           //cycle over chosen frequencies
		doub t,
             maxy=fact*sftab[kk][1],                           //get a size in picture plane for each frequency
			 xg=-maxy+2.*maxy/nxy*doub(iix),                   //offset of current geodesic in x direction
			 yg=-maxy+2.*maxy/nxy*doub(iiy),                   //offset of current geodesic in y direction
			 b=sqrt(xg*xg+yg*yg),                              //impact parameter
			 beta=atan2(yg,xg);                                //polar angle in picture plane


        /*************************************
		 * GEODESIC INTEGRATION FOR EACH RAY *
         *************************************/

		#include "geoint.cpp"
		ppy[currth].nu=1e9*sftab[kk][0];                       //frequency


        /*************************************************
		 * new radiative transfer is solved for each ray *
         *************************************************/

        //clock_t t_b4_solvetrans = clock();
        #include "solvetrans.cpp"
        //t_solvetrans += (clock() - t_b4_solvetrans) / (float)CLOCKS_PER_SEC;


		(*ausin)[iix][iiy][kk][0]=II[0];                       //total intensity I
		(*ausin)[iix][iiy][kk][1]=II[1]*cos(II[2])*sin(II[3]); //linearly polarized intensity Q
		(*ausin)[iix][iiy][kk][2]=II[1]*sin(II[2])*sin(II[3]); //linearly polarized intensity U
		(*ausin)[iix][iiy][kk][3]=II[1]*cos(II[3]);            //circularly polarized intensity V
        //DEFAULT
		(*ausin)[iix][iiy][kk][4]=II[4];                       //any other quantity integrated over a geodesic, i.e., approximate total intensity
		//(*ausin)[iix][iiy][kk][4]=rho; //rho: unknown to this scope...      // column density (rho integrated along a geodesic)

	}; // for(kk=kmin;kk<=kmax;kk+=kstep){

}; // for(ix=0;ix<(nxy+1)*(nxy+1);ix++){
          



//RG:Why "double" and not "doub"?
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
      in[kk][il]*=maxy*maxy;                                              // normalization over integration region size
	totin[kk]=Jy2cgs*ang_size_norm*in[kk][0];                             // normalization for angular size see [global_variables.cpp]
	LPo[kk]=100.*sqrt(in[kk][1]*in[kk][1]+in[kk][2]*in[kk][2])/in[kk][0]; // total LP fraction
	CP[kk]=100.*in[kk][3]/in[kk][0];                                      // CP fraction
	EVPA[kk]=fmod(180/PI*atan2(in[kk][2],in[kk][1])/2.+180.,180);         // EVPA
	err[kk]=Jy2cgs*ang_size_norm*in[kk][4];                               //normalize 5-th intensity as total intensity
	printf("%d; f=%.1f; I=%.3fJy LP=%.2f%% EVPA=%.1fdeg CP=%.3f%% non-pol I=%.2fJy \n",fnum,sftab[kk][0], totin[kk],LPo[kk],EVPA[kk], CP[kk],err[kk]);


    /*****************************
     * WRITE DIAGNOSTICS TO FILE *
     *****************************/

	// setting parameters to be written to file header
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
	sstr<<(int)100*a<<"th"<<(int)floor(100*th+1e-6)<<"f"<<floor(sftab[kk][0])<<"fn"<<fnum<<"case"<<cas<<"_"<<nxy;
    if (SCAN_THROUGH=="view") {
      fif="view";
      sstr<<"_"<<dphi;
    }
	string stra = sstr.str();

	ofstream image_output_file ((dir+"shotimag"+stra+fif+".dat").c_str(), ios::out|ios::binary);   //write into binary output file
	image_output_file.write(reinterpret_cast<char *>(params), 20*sizeof(double));              //write params - fixed size variable
	image_output_file.write(reinterpret_cast<char *>(intab), 5*(nxy+1)*(nxy+1)*sizeof(doub));  //write image
	image_output_file.close();
}

if(iswrite){
    // RG: TODO write header using static int
	stringstream sstr;                                                             //prepare to write full spectra into "poliresa" files
	sstr <<(int)100*a<<"th"<<(int)floor(100*th+1e-6)<<"fn"<<fnum<<"hi";            //"hi": high resolution
	string stra = sstr.str();
	FILE * pFile; 
	pFile=fopen ((dir+"poliresa"+stra+fif+".dat").c_str(),"a");                    //open polires*
	for(kk=kmin;kk<=kmax;kk+=kstep)                                                //write full spectra into "poliresa" files
		fprintf(pFile,"%d %.2f %.5f %.4f %.4f %.4f %.4f %.5f %.4f %.4f %.2f %.4e %.4f\n", fnum,sftab[kk][0],totin[kk],LPo[kk],EVPA[kk],CP[kk],err[kk],heat,rhonor,0.,TpTe,rate*year/Msun,th);
	fclose(pFile);

};

delete [] intab;
delete [] params;
