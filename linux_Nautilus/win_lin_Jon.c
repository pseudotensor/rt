const string dir="/nics/d/home/shcher/analysis/", adir="/lustre/medusa/jmckinne/data1/jmckinne/jmckinne/", fieldstr="/dumps/", xstr="/";
const int nthreads=16;//Nautilus
const int ncuttab[1] = {143}, rlen=272, thlen=128,phlen=256,
	 usgsize=69/*record size in usgdump2d*/,usgoff=156/*offset in usgdump2d file*/;//there is no actual usgdump2d file for these simulations, but the line is needed for consistency with other simulations
const string astr[1]={"thickdisk7"};
const doub atab[1]={0.9375};
