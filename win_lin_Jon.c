// const string dir="/afs/slac.stanford.edu/u/ki/shcher/analysis/", adir="/lustre/ki/pfs/jmckinne/", fieldstr="/dumps/", xstr="/";

//dir  - points to emissivities & other auxilliary files
//adir - points to fluid simulation dump files
//POWERLAW
//const string dir="/home/rgold/rt/powerlaw/", adir="/home/rgold/rt/powerlaw/", fieldstr="/dumps/", xstr="/";

const int nthreads=6; // 24: bh01;// 8;//Orange
const string ASTRORAY_PATH = "/home/rgold/rt/";


// AVERY's TOY+RIAF model (Broderick & Loeb 2009)
// STILL READ IN A MODEL (here a0mad) (AND USE THAT GRID SETUP)
/* const int ncuttab[1] = {120}, rlen=288, thlen=128,phlen=32, */
/*      // actually no usgdump2d file but for consistency */
/* 	 usgsize=69 /\*record size*\/,usgoff=156 /\*offset*\/; */
/* const string astr[1]={""}; */
/* const int dtimdf=5; // stationary model: of no relevance */
/* const doub atab[1]={0.0}; */
/* // a0mad data symbolically linked to toyjet/ */
/* const string dir="/home/rgold/rt/toyjet/", adir=dir, fieldstr="/dumps/", xstr="/"; */



// RG: this is from win_lin_ADAF.c
/* const int ncuttab[3] = {100, 116, 128}, rlen=272, thlen=128,phlen=256, */
/* 	 usgsize=69/\*record size in usgdump2d*\/,usgoff=156/\*offset in usgdump2d file*\/;//there is no actual usgdump2d file for these simulations, but the line is needed for consistency with other simulations */


/******************************
 * RG: THICKDISK7 GRMHD model *
 ******************************/

const string dir=ASTRORAY_PATH+"thermal/", adir=dir, fieldstr="/dumps/", xstr="/";
//RG: WHY "int"?
const int dtimdf=4; // time difference (in M) between fieldlineXXXX files
const int ncuttab[1] = {143}, rlen=272, thlen=128,phlen=256,usgsize=69,usgoff=156;// no actual usgdump2d file,just consistency
const string astr[1]={"thermal"};
const doub atab[1]={0.9375};


// THICKDISKHR3
/* const string dir="/home/rgold/rt/thickdiskhr3/", adir=dir, fieldstr="/dumps/", xstr="/"; */
/* const int dtimdf=4;//time difference (M) between fieldlineXXXX files within consecutive numbers XXXX */
/* const int ncuttab[1] = {143}, rlen=256, thlen=128,phlen=256, */
/*   usgsize=69/\*record size in usgdump2d*\/,usgoff=156; */
/* const string astr[1]={""}; */
/* const doub atab[1]={0.9375}; */


// RG: dipole3dfiduciala GRMHD model
/* const string dir="/home/rgold/rt/dipole3dfiduciala/", adir="/home/rgold/rt/dipole3dfiduciala/", fieldstr="/dumps/", xstr="/"; */

/* const int ncuttab[1] = {120}, rlen=256, thlen=128,phlen=32, */
/* 	 usgsize=69/\*record size in usgdump2d*\/,usgoff=156/\*offset in usgdump2d file*\/;//there is no actual usgdump2d file for these simulations, but the line is needed for consistency with other simulations */
/* const string astr[1]={""}; //"dipole3dfiduciala"}; //RG: POWERLAW */
/* const doub atab[1]={0.92}; */
/* const int dtimdf=2.; */


// RG: quadrupole MRI GRMHD models (i.e. blandford3d)
// http://adsabs.harvard.edu/abs/2009MNRAS.394L.126M
/* const int ncuttab[1] = {120}, rlen=128, thlen=128,phlen=32, */
/*      // actually no usgdump2d file but for consistency */
/* 	 usgsize=69 /\*record size*\/,usgoff=156 /\*offset*\/; */
/* const string astr[1]={""}; */
/* const int dtimdf=4; */
/* const doub atab[1]={0.9375}; */
/* const string dir="/home/rgold/rt/quadrupole/", adir=dir, fieldstr="/dumps/", xstr="/"; */




// RG: thinner MAD GRMHD models
// ORIGINAL
/* const int ncuttab[1] = {120}, rlen=288, thlen=128,phlen=32, */
/* // const int ncuttab[1] = {130}, rlen=288, thlen=128,phlen=32, */
/* // const int ncuttab[1] = {288}, rlen=288, thlen=128,phlen=32, */
/*      // actually no usgdump2d file but for consistency */
/* 	 usgsize=69 /\*record size*\/,usgoff=156 /\*offset*\/; */
/* const string astr[1]={""}; */
/* const int dtimdf=5; */

/* // a0mad or thinnermad0 */
/* const doub atab[1]={0.0}; */
/* const string dir="/home/rgold/rt/a0mad/", adir=dir, fieldstr="/dumps/", xstr="/"; */

// thinnermad9beta25 
/* const doub atab[1]={0.9}; */
/* const string dir="/home/rgold/rt/thinnermad9beta25/", adir=dir, fieldstr="/dumps/", xstr="/"; */

// thinnermadneg9
/* const doub atab[1]={-0.9}; */
/* const string dir="/home/rgold/rt/thinnermadneg9/", adir=dir, fieldstr="/dumps/", xstr="/"; */




/**************************** CHOOSE MODEL **********************************/
/****************************************************************************/

char * descr = getenv("LSB_JOBINDEX");
