const string dir="E:\\GRMHD\\new\\", adir="E:\\GRMHD\\new\\a", fieldstr="\\fieldlines\\", xstr="\\";
const int nthreads=24;//desktop
const int dtimdf=2;//change of time between fieldlineXXXX files within consecutive numbers XXXX
#include "Windows.h"
#include <tchar.h>
const int ncuttab[5] = {100, 110, 116, 128, 138}, rlen=256, thlen=128,phlen=64,
	usgsize=69/*record size in usgdump2d*/,usgoff=156/*offset in usgdump2d file*/;
const string astr[5]={"0","05","07","09","098"};