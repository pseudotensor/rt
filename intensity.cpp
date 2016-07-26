for(ix=0;ix<(snxy+1)*(snxy+1);ix++){                           //cycle over all geodesics (more efficient with OpenMP than cycling along each direction in a double for loop)

    //RG:
    // printf(YELLOW"[intensity.cpp]: "RESET"YO0 ix=%d\n",ix);

	int kk,
		iiy=ix %(snxy+1),                                      //find geodesic coordinates along x and y directions
		iix=(int)ix/(snxy+1);
	for(kk=kmin;kk<=kmax;kk++){                                //cycle over all frequencies
		doub t,
			maxy=fact*sftab[kk][1],                            //get a size in picture plane for each frequency
			xg=-maxy+2.*maxy/snxy*doub(iix),                   //offset of current geodesic in x direction
			yg=-maxy+2.*maxy/snxy*doub(iiy),                   //offset of current geodesic in y direction
			b=sqrt(xg*xg+yg*yg),                               //impact parameter
			beta=atan2(yg,xg);                                 //polar angle in picture plane
		//new geodesic integration for each ray
        
        // doub t_b4_geodesics = clock()/ (doub)CLOCKS_PER_SEC;
        // time_t t_b4_geodesics = time(0);

        //RG:Now moved inside geoint.cpp ~> more general
        //time_t t_b4_geodesics = clock();

        //struct rusage t_b4_geodesics;
        //getrusage(RUSAGE_SELF, &t_b4_geodesics);
        //printf("[intensity.cpp]: t_b4_geodesics=%ld\n",t_b4_geodesics.ru_stime.tv_usec);

        #include "geoint.cpp" // MUST DO FOR EACH FREQUENCY ONLY BECAUSE IMAGE PLANE SIZE IS FREQUENCY DEPENDENT

        //RG:Now moved inside geoint.cpp ~> more general
        // t_geodesics += (clock() - t_b4_geodesics) / (doub)CLOCKS_PER_SEC;

        // time_t t_after_geodesics = time(0);
        // time_t t_after_geodesics = time(0);
        // t_geodesics += difftime(t_after_geodesics,t_b4_geodesics);
        //getrusage(RUSAGE_SELF, &t_after_geodesics);
        //t_geodesics += (t_after_geodesics.ru_stime.tv_usec - t_b4_geodesics.ru_stime.tv_usec);

		ppy[currth].nu=1e9*sftab[kk][0];                       //frequency

        // printf(YELLOW"[intensity.cpp]: "RED"PRE [solvetrans.cpp]"RESET"\n");
        //doub t_b4_solvetrans = clock();
        //clock_t t_b4_solvetrans = clock();
		//new radiative transfer is solved for each ray
        #include "solvetrans.cpp"
        // t_solvetrans += (clock() - t_b4_solvetrans) / (float)CLOCKS_PER_SEC;
        // printf(YELLOW"[intensity.cpp]: "RED"POST [solvetrans.cpp]"RESET"\n");
		
		(*ausin)[iix][iiy][kk][0]=II[0];                       //total intensity I
		(*ausin)[iix][iiy][kk][1]=II[1]*cos(II[2])*sin(II[3]); //linearly polarized intensity Q
		(*ausin)[iix][iiy][kk][2]=II[1]*sin(II[2])*sin(II[3]); //linearly polarized intensity U
		(*ausin)[iix][iiy][kk][3]=II[1]*cos(II[3]);            //circularly polarized intensity V
        for (int idx=4; idx<nr_of_imagediags;idx++)
          (*ausin)[iix][iiy][kk][idx]=II[idx];                       //any other quantity integrated over a geodesic, i.e., approximate total intensity
	};
};

doub in[sflen][5];                                             //define array for computing intensity

for(kk=kmin;kk<=kmax;kk++){                                    

	for(il=0;il<5;il++) in[kk][il]=0.;                         // INITIALIZE TO ZERO

	doub hei=2./snxy;                                          //distance between neighbooring points in picture plane //RG: image plane side length = 2?
    doub maxy=fact*sftab[kk][1];                               //size of the integration regions in picture plane 

    if (INTERPOLATION_METHOD=="SUM") {
      for(ix=0;ix<=snxy;ix+=1)                                 //2D integration loop over picture plane
        for(iy=0;iy<=snxy;iy+=1)
          for(il=0;il<5;il++) {                               //for each kind of intensity [I,Q,U,V,5th diagnostic]

            //RG: CHECK BY TRYING SEVERAL FROM THIS LINK AND COMPARE
            // https://www.kth.se/social/upload/52a04c17f27654620e188bb0/Ant-Integration.pdf

            if ((il+ix+iy==0) && (kk==kmin)) printf(YELLOW"[intensity.cpp]: "RESET"INTERPOLATION METHOD=SUM\n");
            in[kk][il]+=(*ausin)[ix][iy][kk][il]*hei*hei; // ASSUME CONST VALUE WITHIN CELLS (0th ORDER, i.e. O(dx) APPROXIMATION)
          }
    }

    else if (INTERPOLATION_METHOD=="ASTRORAYv1.0") {

      for(ix=0;ix<=snxy-3;ix+=2)                              //2D integration loop over picture plane
        for(iy=0;iy<=snxy-3;iy+=2)
          for(il=0;il<5;il++) {                               //for each kind of intensity [I,Q,U,V,5th diagnostic]


            // 2nd order 2D stencil shifted:
            //    ~> j
            // |   / 1  4 1 \ 
            // v  (  4 16 4  )
            // i   \ 1  4 1 / 

            // 2nd order accuracy of the integrator - works best RG: What does best mean? In what way?
            in[kk][il]+=((*ausin)[ix][iy][kk][il] + (*ausin)[ix][2 + iy][kk][il] + 4*(*ausin)[1 +ix][iy][kk][il] + 4*(*ausin)[1 + ix][2 + iy][kk][il] + (*ausin)[2 +ix][iy][kk][il] +
                         4*((*ausin)[ix][1 + iy][kk][il] + 4*(*ausin)[1 + ix][1 +iy][kk][il] + (*ausin)[2 + ix][1 + iy][kk][il]) + (*ausin)[2 + ix][2 +iy][kk][il])*hei*hei/9.;
          } // for x for y for il


            // BOUNDARY

            for(ix=0;ix<=snxy-1;ix++)                                  //integration over x boundary - 1-st order
              for(il=0;il<5;il++) {
                in[kk][il]+=((*ausin)[ix][snxy-1][kk][il] + (*ausin)[ix][snxy][kk][il] + (*ausin)[1 + ix][snxy-1][kk][il] + (*ausin)[1 + ix][snxy][kk][il])*hei*hei/4.;
              }
            for(iy=0;iy<=snxy-2;iy++)                                  //integration over y boundary - 1-st order
              for(il=0;il<5;il++) {
                in[kk][il]+=((*ausin)[snxy-1][iy][kk][il] + (*ausin)[snxy-1][1 + iy][kk][il] + (*ausin)[snxy][iy][kk][il] + (*ausin)[snxy][1 + iy][kk][il])*hei*hei/4.;
              }
    } // else if (INTERPOLATION_METHOD=="ASTRORAYv1.0") {

	for(il=0;il<5;il++) in[kk][il]*=maxy*maxy;                            // normalization over integration region size  

	totin[kk]=Jy2cgs*ang_size_norm*in[kk][0];                             // normalization for angular size of Sgr A* region - recompute for your object!
	LPo[kk]=100.*sqrt(in[kk][1]*in[kk][1]+in[kk][2]*in[kk][2])/in[kk][0]; // compute total LP fraction <LP> = <Q,U>/<I>
	CP[kk]=100.*in[kk][3]/in[kk][0];                                      // compute CP fraction       <CP> = <V>  /<I>
	EVPA[kk]=fmod(180/PI*atan2(in[kk][2],in[kk][1])/2.+180.,180.);        // compute EVPA
	err[kk]=Jy2cgs*ang_size_norm*in[kk][4];                               // normalize 5-th intensity as total intensity
	printf(YELLOW"[intensity.cpp]:"RESET" fnum=%d; f=%.1f; I=%.3fJy LP=%.2f%% EVPA=%.1fdeg CP=%.3f%% non-pol I=%.2fJy \n",fnum,sftab[kk][0], totin[kk],LPo[kk],EVPA[kk], CP[kk],err[kk]);
};


// Chi^2
//RG:FLAG Chi^2/dof should be computed by a single fct! Different (w/o nu[5-10]) Chi^2 in m_search e.g.
// ~> FCT WRITTEN 2B INCORPORATED... Oct 21st 2015
doub chi=0.,chi_I=0.;doub resid[14]; //[4];
chisquare(totin,LPo,CP,chi,chi_I,resid); // ~> see [chisquare.cpp]
printf(YELLOW"[intensity.cpp]: "MAGENTA"chi=%f,chi_I=%f"RESET"\n",chi,chi_I);

// Hm something is different...
// chisquare(): chi^2/dof=9.856375 chi^2_I/dof=13.966356
//              [intensity.cpp]: chi=9.856375,chi_I=13.966356
// # Chi^2/dof = 8.267417 [Based on SED only] poliresa file

doub xisq,                                                               //\chi^2
	 dof=7.;                                                             //degrees of freedom //RG:FLAG:HARDCODED

// totin[4] ~> F@84Ghz
xisq=pow((doub)(totin[4]-tofit[4][1])/dFnu[4],2)+pow((doub)(totin[5]-tofit[5][1])/dFnu[5],2)+pow((doub)(totin[6]-tofit[6][1])/dFnu[6],2)+
	 pow((doub)(totin[7]-tofit[7][1])/dFnu[7],2)+pow((doub)(totin[8]-tofit[8][1])/dFnu[8],2)+pow((doub)(totin[9]-tofit[9][1])/dFnu[9],2)+
	 pow((doub)(totin[10]-tofit[10][1])/dFnu[10],2);
xisq/=dof;//compute reduced \chi^2

if(iswrite){
	stringstream sstr;                                                  //prepare to write full spectra into "poliresa" files
	// sstr<<(int)100*a<<"th"<<(int)floor(100*th+1e-6)<<"fn"<<fnum<<"boo"; //"boo" is a legacy string
	sstr<<(int)100*a<<"th"<<(int)floor(100*th+1e-6)<<"fn"<<fnum<<"case"<<cas<<"boo"; //"boo" is a legacy string
	string stra = sstr.str();
	FILE * pFile;
	pFile=fopen ((dir+"poliresa"+stra+fif+".dat").c_str(),"a");         //"fopen" is unsafe, but it works. If done w/ streams, then no flexibility of "fprintf"
	if(pFile==NULL){
		printf("Cannot open poliresa output file \n Exiting \n");
		exit(-1);
	};

    // HEADER (with xisq)
    fprintf(pFile,"# Chi^2/dof = %f [Based on SED only]\n# fnum,sftab[kk][0],totin[kk],LPo[kk],EVPA[kk],CP[kk],err[kk],heat,rhonor,xisq,TpTe,rate*year/Msun,th \n",xisq);

	for(kk=kmin;kk<=kmax;kk++)                                          //write full spectra into "poliresa" files
		fprintf(pFile,"%d %.2f %.5f %.3f %.3f %.3f %.3f %.5f %.4f %.4f %.2f %.4e %.4f\n",fnum,sftab[kk][0],totin[kk],LPo[kk],EVPA[kk],CP[kk],err[kk],heat,rhonor,xisq,TpTe,rate*year/Msun,th);
	fclose(pFile);
};
