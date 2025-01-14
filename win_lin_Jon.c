const string dir="/afs/slac.stanford.edu/u/ki/shcher/analysis/", adir="/lustre/ki/pfs/jmckinne/", fieldstr="/dumps/", xstr="/";
const int nthreads=8;//Orange
const int dtimdf=4;//time difference (M) between fieldlineXXXX files within consecutive numbers XXXX
const int ncuttab[1] = {143}, rlen=272, thlen=128,phlen=256,
	 usgsize=69/*record size in usgdump2d*/,usgoff=156/*offset in usgdump2d file*/;//there is no actual usgdump2d file for these simulations, but the line is needed for consistency with other simulations
const string astr[1]={"thickdisk7"};
const doub atab[1]={0.9375};
char * descr = getenv("LSB_JOBINDEX");
