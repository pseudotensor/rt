{ 
  // obtain averages of temperature and density profiles to find Te(Ts) and Tp(Ts) functions
  // This routine only averages simulation data over a specified range in snapshots
  // It therefore only needs to be rerun when any of the following changes: simulation, radial extension (dxdxp.dat), snapshot interval (fmin,fmax)
  int 
    nx1,nx2,nx3,                  //grid size along r, theta, and phi directions in fluid simulation snapshots
	fmin=atoi(argv[2]),              //minimum ID of fluid simulation snapshot = 2nd command line argument
    fmax=atoi(argv[3]),             //maximum ID of fluid simulation snapshot = 3nd command line argument
  //fmax=atoi(),             //maximum ID of fluid simulation snapshot = 3nd command line argument
	sep=1,                //ID difference between consecutive considered fluid simulation snapshots = 1 => each snapshot is considered
	ind=(fmax-fmin)/sep+1;//total number of snapshots
  doub 
    dx1, off, x1min, tt, //auxiliary quantites read as text (not as binary) from fluid simulation file - check your fluid simulation file format!
    rad[rlen],           //radial grid read from auxiliary "usgdump2d" or "dxdxp.dat" file - check your fluid simulation file format!
    rho[rlen],           //variable for computing mean (equatorial) radial density profile
    ts[rlen],            //variable for computing mean (equatorial) radial energy density profile
    a=atab[sp], asq=a*a; // define spin value a and a^2

typedef double (*usgarr)[thlen][rlen][usgsize];//set of records for each point on 2D fluid simulation grid
usgarr usgread = (usgarr) new double[thlen][rlen][usgsize];

// Following variables only initialized to avoid problems with uninitialized variables
rhonor=1000000.; heat=0.5; th=1.0; fdiff=0;



/************************************************
 * GET SIMULATION PARAMETERS: time & resolution *
 ************************************************/

string fdir=adir+astr[sp]+fieldstr; // directory name for fluid simulation snapshot files
stringstream sstr;                
sstr<<setfill('0')<<setw(4)<<fmin;                                    //build filename of a specific fluid simulation snapshot
cout << YELLOW"[m_ts.cpp]: "RESET"Assuming fluid snapshot is here: fdir="+fdir+"fieldline"+sstr.str()+".bin"<<endl;
ifstream xfline((fdir+"fieldline"+sstr.str()+".bin").c_str(),ios::in);//read from a fluid simulation snapshot as text
// USE NEW PREPROCESSOR DIRECTIVES IF MODELS DIFFER
xfline>>tt>>nx1>>nx2>>nx3>>x1min>>dx1>>dx1>>dx1;                      //read a number of variables from file including grid dimensions - check how it works for your fluid simulations
xfline.close();

printf(YELLOW"[m_ts.cpp]: "RESET"tt=%e,nx1=%d,nx2=%d,nx3=%d,x1min=%e,dx1=%e\n",tt,nx1,nx2,nx3,x1min,dx1); // talk about it

if((rlen!=nx1)||(thlen!=nx2)||(phlen!=nx3)) { // check that the grid dimensions specified in "win_lin..." agrees with that found in simulation snapshots
  printf(YELLOW"[m_ts.cpp]: "RESET"rlen=%d,nx1=%d,thlen=%d,nx2=%d,phlen=%d,nx3=%d\nErrors in dimensions \n Exiting ",rlen,nx1,thlen,nx2,phlen,nx3);
  exit(-1);
};



/******************************************
 * GET COORDINATES INCL. RADIAL EXTENSION *
 ******************************************/

filebuf *pbuf;                                                        //auxiliary file buffer
//if(astr[sp].length()<4){                                              //differentiate between reading "usgdump2d" and "dxdxp.dat"

 // if(false){                                              //differentiate between reading "usgdump2d" and "dxdxp.dat"
 if(RADIAL_EXT_FILE=="usgdump2d"){                                              //differentiate between reading "usgdump2d" and "dxdxp.dat"
    cout << "Assuming usgdump2d.dat is here: "+dir+astr[sp]+xstr<<endl;
	ifstream gre((dir+astr[sp]+xstr+"usgdump2d").c_str(), ios::in|ios::binary);//read "usgdump2d" file as binary - check path!
	pbuf=gre.rdbuf();                                                 //define buffer
	pbuf->pubseekpos(usgoff);                                         //set the reading position at the beginning of the array of records
	int tosize=usgsize*rlen*thlen*sizeof(double);                     //compute the size of the array of records
	gre.read(reinterpret_cast<char *>(*usgread), tosize);             //read the array of records
	gre.close();
	for(int k=0;k<rlen;k++){
		rad[k]=(*usgread)[0][k][7];                                   //define radial grid (independent of theta)
	};
	for(int k=0;k<rlen;k++)
		for(int i=0;i<thlen;i++)
			theta[k][i]=-cos((*usgread)[i][k][8]);                    //define theta grid (for each radius)
} else {
    cout << "Assuming dxdxp.dat is here: "+dir+astr[sp]+xstr<<endl;
	ifstream dxp((dir+astr[sp]+xstr+"dxdxp.dat").c_str(), ios::in|ios::binary); // read "dxdxp.dat" file, which was pre-generated in Mathematica
	pbuf=dxp.rdbuf();                                                           // define buffer
	dxp.read(reinterpret_cast<char *>(coord), ndd*thlen*2*sizeof(float));       // read coordinates 2D matrix
	dxp.read(reinterpret_cast<char *>(dxdxp), ndd*thlen*4*4*sizeof(float));     // read transformation matrix from MKS to KS - not actually used // RG: WHAT?!
	dxp.close();

	for(int k=0;k<rlen;k++)
      for(int i=0;i<thlen;i++){
        if (isnan(coord[k][i][1])) {
          cout << "[m_ts.cpp]: Theta coordinates contain nan. Unwise to continue..." << endl;
          exit(1);
        }
        theta[k][i]=coord[k][i][1];                               //define theta grid (for each radius)
      }

	for(int k=0;k<rlen;k++){
      if (isnan(coord[k][0][0])) {
        cout << "[m_ts.cpp]: Radial coordinates contain nan. Unwise to continue..." << endl;
        exit(1);
      }
      rad[k]=exp(coord[k][0][0]);                                   //radial grid
	};
 }



/********************
 * COMPUTE AVERAGES *
 ********************/

// initialize density and energy density radial profiles
for(int rk=0;rk<rlen;rk++){
  rho[rk]=0.;
  ts[rk]=0.;
}
doub rhoinst=0.;
int nx=37; // RG:HARDWIRED AS IN ASTRORAYv1.0; BAD BUT JUST A ROUGH DIAGNOSTIC... // radial index, where instantaneous density-average is computed

//RG:TODO can openMP this, BUT need to make rho,ts,uu thread-private
for(fnum=fmin;fnum<=fmax;fnum+=sep){                                  //cycle over fluid simulation snapshots
	stringstream sstr;
	sstr<<setfill('0')<<setw(4)<<fnum;
	ifstream fline((fdir+"fieldline"+sstr.str()+".bin").c_str(),ios::in|ios::binary);//read another fluid simulation snapshots
	pbuf=fline.rdbuf();                                               //initialize buffer
	int fsize=pbuf->pubseekoff (0,ios::end),                          //size of file
	tosize=wdd*phlen*thlen*rlen*sizeof(float);                        //size of binary section of fluid simulation snapshot file
	pbuf->pubseekpos(fsize-tosize);                                   //set file position at the beginning of binary section
	fline.read(reinterpret_cast<char *>(uu[0]), tosize);              //read fluid simulation snapshot
	fline.close();

    //RG:TODO can openMP this trivially
	for(int rk=0;rk<rlen;rk++) {
      for(int thk=thlen/2-3;thk<=thlen/2+3;thk++) {                   // only over theta angles close to the equatorial plane
        for(int phk=0;phk<phlen;phk++) {
              rho[rk]+=(*uu[0])[phk][thk][rk][0]/7/phlen/ind; // get average over phi and +/-3 cells around equator of density (rho) at each radius
              ts[rk]+=(*uu[0])[phk][thk][rk][1]/7/phlen/ind*mp*cc*cc/3/kb/(*uu[0])[phk][thk][rk][0]; // get mean energy density (u/rho) at each radius
              if (rk==nx) // RG: HARDWIRED AS IN ASTRORAYv1.0... NOT GOOD BUT JUST A ROUGH DIAGNOSTIC ANYWAY...
                rhoinst+=(*uu[0])[phk][thk][nx-1][0]/7/phlen;             //instantaneous space-averaged density
        }
      }            
      doub rx=rad[nx-1];                                                // radius corresponding to nx
      // printf(YELLOW"[m_ts.cpp]: "RED"[HARDWIRE]: nx=%d rx=%f"RESET"\n",nx,rx);
      if (rk==nx) printf(YELLOW"[m_ts.cpp]: "RESET"fn=%d: average density at r[%d-1]=%f is rho=%f \n",fnum,nx,rx,rhoinst);
    }
}



/************************
 * BUBBLE SORTING ARRAY *
 ************************/

bool sorted=false;                                                    //bubble sorting of density and energy density profiles. These arrays are already almost sorted
doub temp;
while(!sorted){
    sorted=true;
	for(ix=1;ix<rlen;ix++){
		if(ts[ix-1]<ts[ix]){
			temp=ts[ix-1];
			ts[ix-1]=ts[ix];
			ts[ix]=temp;
			sorted=false;
		};
		if(rho[ix-1]<rho[ix]){
			temp=rho[ix-1];
			rho[ix-1]=rho[ix];
			rho[ix]=temp;
			sorted=false;
		};
	};
};


/******************
 * OUTPUT TO FILE *
 ******************/

FILE * yfline;                                                       //file variable
yfline=fopen((adir+astr[sp]+xstr+"Tsmapx"+astr[sp]+".dat").c_str(),"w");//Tsmapx file
for(int rk=0;rk<rlen;rk++)
	fprintf(yfline,"%.8e %.8e %.8e \n",rad[rk],ts[rk],rho[rk]);      //writing computed radial profiles of density/energy density into file
fclose(yfline);
}
