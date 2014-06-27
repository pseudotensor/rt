{//averages temperature and density profiles to find Te(Ts) and Tp(Ts) functions
int i, k, 
	nx1,                  //grid size along r, theta, and phi directions in fluid simulation snapshots
	nx2,
	nx3,
	rk,                   //array indices in r, theta, and phi directions
	thk,
	phk,
	fmin=co,              //minimum ID of fluid simulation snapshot = 2nd command line argument
	fmax=mco,             //maximum ID of fluid simulation snapshot = 3nd command line argument
	sep=1,                //ID difference between consecutive considered fluid simulation snapshots = 1 => each snapshot is considered
	ind=(fmax-fmin)/sep+1;//total number of snapshots
doub dx1, off, x1min, tt, //auxiliary quantites read as text (not as binary) from fluid simulation file - check your fluid simulation file format!
	 rad[rlen],           //radial grid read from auxiliary "usgdump2d" file - check your fluid simulation file format!
	 rho[rlen],           //variable for computing mean (equatorial) radial density profile
	 ts[rlen];            //variable for computing mean (equatorial) radial energy density profile
typedef double (*usgarr)[thlen][rlen][usgsize];//set of records for each point on 2D fluid simulation grid
usgarr usgread = (usgarr) new double[thlen][rlen][usgsize];

rhonor=1000000.;          //initialize a sample model to avoid problems with uninitialized variables
heat=0.5;
th=1.0;
fdiff=0;
                          
string fdir=adir+astr[sp]+fieldstr;//directory name for fluid simulation snapshot files
a=atab[sp];                        //define spin value a and a^2
asq=a*a;
stringstream sstr;                
sstr<<setfill('0')<<setw(4)<<fmin;                                    //build filename of a specific fluid simulation snapshot
ifstream xfline((fdir+"fieldline"+sstr.str()+".bin").c_str(),ios::in);//read from a fluid simulation snapshot as text
xfline>>tt>>nx1>>nx2>>nx3>>x1min>>dx1>>dx1>>dx1;                      //read a number of variables from file including grid dimensions - check how it works for your fluid simulations
xfline>>off>>off>>off>>off>>off>>off;
xfline>>off>>off>>off;
xfline.close();
if((rlen!=nx1)||(thlen!=nx2)||(phlen!=nx3)){                          //check that the grid dimensions specified in "win_lin..." agrees with that found in simulation snapshots
	printf("Errors in dimensions \n Exiting ");
	exit(-1);
};


filebuf *pbuf;                                                        //auxiliary file buffer
if(astr[sp].length()<4){                                              //differentiate between reading "usgdump2d" and "dxdxp.dat"
	ifstream gre((dir+astr[sp]+xstr+"usgdump2d").c_str(), ios::in|ios::binary);//read "usgdump2d" file as binary - check path!
	pbuf=gre.rdbuf();                                                 //define buffer
	pbuf->pubseekpos(usgoff);                                         //set the reading position at the beginning of the array of records
	int tosize=usgsize*rlen*thlen*sizeof(double);                     //compute the size of the array of records
	gre.read(reinterpret_cast<char *>(*usgread), tosize);             //read the array of records
	gre.close();
	for(i=0;i<rlen;i++){
		rad[i]=(*usgread)[0][i][7];                                   //define radial grid (independent of theta)
	};
	for(k=0;k<rlen;k++)
		for(i=0;i<thlen;i++)
			theta[k][i]=-cos((*usgread)[i][k][8]);                    //define theta grid (for each radius)
} else {
	ifstream dxp((dir+astr[sp]+xstr+"dxdxp.dat").c_str(), ios::in|ios::binary);//read "dxdxp.dat" file, which was pre-generated in Mathematica
	pbuf=dxp.rdbuf();                                                 //define buffer
	dxp.read(reinterpret_cast<char *>(coord), ndd*thlen*2*sizeof(float));//read coordinates 2D matrix
	dxp.read(reinterpret_cast<char *>(dxdxp), ndd*thlen*4*4*sizeof(float));//read transformation matrix from MKS to KS - not actually used
	dxp.close();
	for(k=0;k<rlen;k++)
		for(i=0;i<thlen;i++){
			theta[k][i]=coord[k][i][1];                               //define theta grid (for each radius)
		}
	for(k=0;k<rlen;k++){
		rad[k]=exp(coord[k][0][0]);                                   //radial grid
	};
}

for(rk=0;rk<rlen;rk++){                                               //initialize density and energy density radial profiles
	rho[rk]=0.;
	ts[rk]=0.;
}

for(fnum=fmin;fnum<=fmax;fnum+=sep){                                  //cycle over fluid simulation snapshots
	stringstream sstr;
	sstr<<setfill('0')<<setw(4)<<fnum;
	ifstream fline((fdir+"fieldline"+sstr.str()+".bin").c_str(),ios::in|ios::binary);//read another fluid simulation snapshots
	pbuf=fline.rdbuf();                                               //initialize buffer
	int fsize=pbuf->pubseekoff (0,ios::end),                          //size of file
		tosize=11*phlen*thlen*rlen*sizeof(float);                     //size of binary section of fluid simulation snapshot file
	pbuf->pubseekpos(fsize-tosize);                                   //set file position at the beginning of binary section
	fline.read(reinterpret_cast<char *>(uu[0]), tosize);              //read fluid simulation snapshot
	fline.close();
	for(rk=0;rk<rlen;rk++)                                            //for all radii
		for(thk=thlen/2-3;thk<=thlen/2+3;thk++)                       //for theta angles close to the equatorial plane
			for(phk=0;phk<phlen;phk++){                               //for all azimuthal angles
				rho[rk]+=(*uu[0])[phk][thk][rk][0]/7/phlen/ind;       //compute mean density at that radius
				ts[rk]+=(*uu[0])[phk][thk][rk][1]/7/phlen/ind*mp*cc*cc/3/kb/(*uu[0])[phk][thk][rk][0];//compute mean energy density at that radius
			};

//done till here
	int nx=37;                                                        //radial point, where instantaneous space-averaged density is computed
	doub rx=rad[nx-1];                                                //correspondent radius
	doub rhoinst=0.;

	for(thk=thlen/2-3;thk<=thlen/2+3;thk++)                           //for theta angles close to the equatorial plane
		for(phk=0;phk<phlen;phk++)                                    //for all azimuthal angles
			rhoinst+=(*uu[0])[phk][thk][nx-1][0]/7/phlen;             //instantaneous space-averaged density
	printf("fn=%d: average density at r=%f is rho=%f \n",fnum,rx,rhoinst);
};
	
bool sorted=false;                                                    //bubble sorting of density and energy density profiles. These arrays are alredy almost sorted
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

FILE * yfline;                                                       //file variable
yfline=fopen((adir+astr[sp]+xstr+"Tsmapx"+astr[sp]+".dat").c_str(),"w");//Tsmapx file
for(rk=0;rk<rlen;rk++)
	fprintf(yfline,"%.8e %.8e %.8e \n",rad[rk],ts[rk],rho[rk]);      //writing computed radial profiles of density/energy density into file
fclose(yfline);
}