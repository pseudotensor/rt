const string dir="E:\\GRMHD\\ADAFs\\", adir="E:\\GRMHD\\ADAFs\\", fieldstr="\\fieldlines\\", xstr="\\";
const int nthreads=24;//desktop
const int dtimdf=2;//change of time between fieldlineXXXX files within consecutive numbers XXXX
#include "Windows.h"
#include <tchar.h>
const int ncuttab[5] = {100, 110, 116, 128, 138}, rlen=264, thlen=126,phlen=60,
	usgsize=69/*record size in usgdump2d*/,usgoff=156/*offset in usgdump2d file*/;
const string astr[5]={"a0","a05","a07","a09","a098"};
const doub atab[5]={0., 0.5, 0.7, 0.9, 0.98};