//simulations by Olek on Odyssey blackhole
const string dir="/afs/slac.stanford.edu/u/ki/shcher/analysis/", adir="/lustre/ki/pfs/jmckinne/", fieldstr="/dumps/", xstr="/";//Orange cluster at Stanford
//const string dir="E:\\GRMHD\\", adir="E:\\GRMHD\\", fieldstr="\\fieldlines\\", xstr="\\";                                      //Windows workstation at UMD
//dir  - points to emissivities & other auxilliary files
//adir - points to fluid simulation dump files
const int nthreads=8;                           //Orange
const int dtimdf=10;                            //time difference between fieldlineXXXX files within consecutive numbers XXXX and XXXX+1
const int ncuttab[3] = {100, 116, 128},         //the last radial grid point, where the simulation is considered converged
	rlen=264, thlen=126, phlen=60,           //grid dimensions of fluid dynamics simulations in polar coordinates
	usgsize=69,                              //record size in usgdump2d file
	usgoff=156;                              //offset in usgdump2d file
const doub atab[3]={0., 0.7, 0.9};              //BH spin values
const string astr[3]={"Oleka0","Oleka07","Oleka09"};//subdirectories for each spin

//MVS thinks getenv is unsafe, but hard to make it platform-independent otherwise
//environment variable PBS_ARRAYID (QSUB scheduler) or LSB_JOBINDEX (BSUB scheduler) is the job ID, when submitted as a job array
//when debugging/running on workstation, one has to explicitly create this environment variable
char * descr = getenv("LSB_JOBINDEX");
