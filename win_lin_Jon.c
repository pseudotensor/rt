// const string dir="/afs/slac.stanford.edu/u/ki/shcher/analysis/", adir="/lustre/ki/pfs/jmckinne/", fieldstr="/dumps/", xstr="/";

//dir  - points to emissivities & other auxilliary files
//adir - points to fluid simulation dump files
//POWERLAW
//const string dir="/home/rgold/rt/powerlaw/", adir="/home/rgold/rt/powerlaw/", fieldstr="/dumps/", xstr="/";
//const string MODEL="dipole"; // dipole,quadrupole,thickdisk7
const int nthreads=24; // 24: bh01;// 8;//Orange
const string ASTRORAY_PATH = "/home/rgold/rt/";
#define THINNERMAD0 0
#define THICKDISK7 7
#define DIPOLE 2
#define QUADRUPOLE 4
#define TOYJET 100

/* USER CHOICE HERE */
#define MODEL DIPOLE


#if MODEL == THICKDISK7
// AVERY's TOY+RIAF model (Broderick & Loeb 2009)
// STILL READ IN A MODEL (here a0mad) (AND USE THAT GRID SETUP)
const int ncuttab[1] = {120}, rlen=288, thlen=128,phlen=32,
     // actually no usgdump2d file but for consistency
	 usgsize=69 /*record size*/,usgoff=156 /*offset*/;
const string astr[1]={""};
const int ndd=350;           //radial dimension of coordinate/coordinate transformation matrices
const int dtimdf=5; // stationary model: of no relevance
const doub atab[1]={0.0};
// a0mad data symbolically linked to toyjet/
const string dir="/home/rgold/rt/toyjet/", adir=dir, fieldstr="/dumps/", xstr="/";

/* const int ncuttab[3] = {100, 116, 128}, rlen=272, thlen=128,phlen=256, */ // win_lin_ADAF.c



/******************************
 * RG: THICKDISK7 GRMHD model *
 ******************************/
#elif MODEL == THICKDISK7
  const string dir=ASTRORAY_PATH+"thermal/",
    // const string dir=ASTRORAY_PATH+"thickdisk7/",
    adir=dir, fieldstr="/dumps/", xstr="/";
  //RG: WHY "int"?
  const int ndd=650;           //radial dimension of coordinate/coordinate transformation matrices
  const int dtimdf=4; // time difference (in M) between fieldlineXXXX files
  const int ncuttab[1] = {143}, rlen=272, thlen=128,phlen=256,usgsize=69,usgoff=156;// no actual usgdump2d file,just consistency
  const string astr[1]={"thermal"};
  // const string astr[1]={"thickdisk7"};
  const doub atab[1]={0.9375}; 

#elif MODEL == THICKDISKHR3
// THICKDISKHR3
/* const string dir="/home/rgold/rt/thickdiskhr3/", adir=dir, fieldstr="/dumps/", xstr="/"; */
/* const int ndd=350;           //radial dimension of coordinate/coordinate transformation matrices */
/* const int dtimdf=4;//time difference (M) between fieldlineXXXX files within consecutive numbers XXXX */
/* const int ncuttab[1] = {143}, rlen=256, thlen=128,phlen=256, */
/*   usgsize=69/\*record size in usgdump2d*\/,usgoff=156; */
/* const string astr[1]={""}; */
/* const doub atab[1]={0.9375}; */

#elif MODEL == DIPOLE
// RG: dipole3dfiduciala GRMHD model

const string dir="/home/rgold/rt/dipole3dfiduciala/", adir="/home/rgold/rt/dipole3dfiduciala/", fieldstr="/dumps/", xstr="/";
const int ncuttab[1] = {120}, rlen=256, thlen=128,phlen=32,
	 usgsize=69/*record size in usgdump2d*/,usgoff=156/*offset in usgdump2d file*/;//there is no actual usgdump2d file for these simulations, but the line is needed for consistency with other simulations
const string astr[1]={""}; //"dipole3dfiduciala"}; //RG: POWERLAW
const doub atab[1]={0.92};
const int ndd=350;           //radial dimension of coordinate/coordinate transformation matrices
const int dtimdf=2.;



#elif MODEL == QUADRUPOLE
  // RG: quadrupole MRI GRMHD models (i.e. blandford3d)
  // http://adsabs.harvard.edu/abs/2009MNRAS.394L.126M
  const int ncuttab[1] = {120}, rlen=128, thlen=128,phlen=32,
     // actually no usgdump2d file but for consistency
    usgsize=69 /*record size*/,usgoff=156 /*offset*/;
  const string astr[1]={""};
  const int ndd=350;           //radial dimension of coordinate/coordinate transformation matrices
  const int dtimdf=2;
  const doub atab[1]={0.9375};
  const string dir="/home/rgold/rt/quadrupole/", adir=dir, fieldstr="/dumps/", xstr="/";



//#elif MODEL == THINNERMAD
#elif MODEL == THINNERMAD0
// RG: thinner MAD GRMHD models
// ORIGINAL
const int ncuttab[1] = {120}, rlen=288, thlen=128,phlen=32,
// const int ncuttab[1] = {130}, rlen=288, thlen=128,phlen=32,
// const int ncuttab[1] = {288}, rlen=288, thlen=128,phlen=32,
     // actually no usgdump2d file but for consistency
	 usgsize=69 /*record size*/,usgoff=156 /*offset*/;
const string astr[1]={""};
const int ndd=350;           //radial dimension of coordinate/coordinate transformation matrices
const int dtimdf=5;

/* // a0mad or thinnermad0 */
const doub atab[1]={0.0};
const string dir="/home/rgold/rt/a0mad/", adir=dir, fieldstr="/dumps/", xstr="/";

// RG: define other models or have separate spin variable?
// thinnermad9beta25 
/* const doub atab[1]={0.9}; */
/* const string dir="/home/rgold/rt/thinnermad9beta25/", adir=dir, fieldstr="/dumps/", xstr="/"; */

// thinnermadneg9
/* const doub atab[1]={-0.9}; */
/* const string dir="/home/rgold/rt/thinnermadneg9/", adir=dir, fieldstr="/dumps/", xstr="/"; */

#endif

char * descr = getenv("LSB_JOBINDEX");


// some experiments/snippets with tenary operator
// int ndd;
// THICKDISK7
// if ( !strcmp(astr[0],"thickdisk7") ) {

// if (true) {
// if (!strcmp(dir,"thickdisk7") {

// const int i = someCondition ? calculatedValue : defaultValue;
// DOES NOT WORK
//const int ndd = (astr[0].c_str()=="thickdisk7") ? 650 : 350;
//const string astr[1]={"thickdisk7"};

// THIS WORKS ON http://www.tutorialspoint.com/compile_cpp_online.php but not here...
// const int ndd= astr[0]=="thickdisk7" ? 650 : 350;

// WORKS
// const int ndd = (true) ? 650 : 350;
