const string dir="E:\\GRMHD\\", adir="E:\\GRMHD\\", fieldstr="\\fieldlines\\", xstr="\\";
const int nthreads=48;//desktop
const int dtimdf=4;//change of time between fieldlineXXXX files within consecutive numbers XXXX
#include "Windows.h"
#include <tchar.h>
const int ncuttab[1] = {143}, rlen=272, thlen=128,phlen=256,
	 usgsize=69/*record size in usgdump2d*/,usgoff=156/*offset in usgdump2d file*/;//there is no actual usgdump2d file for these simulations, but the line is needed for consistency with other simulations
const string astr[1]={"thickdisk7"};
const doub atab[1]={0.9375};
char * descr = getenv("LSB_JOBINDEX");