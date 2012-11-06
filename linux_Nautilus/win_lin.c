const string dir="E:\\GRMHD\\", adir="E:\\GRMHD\\", fieldstr="\\fieldlines\\", xstr="\\";
const int nthreads=24;//desktop
#include "Windows.h"
#include <tchar.h>
const int ncuttab[5] = {145, 147, 149, 153, 162}, rlen=256, thlen=64,phlen=32, 
	usgsize=70/*record size in usgdump2d*/,usgoff=168/*offset in usgdump2d file*/;
const string astr[5]={"a0","a05","a07","a09","a098"};
const doub atab[5]={0., 0.5, 0.7, 0.9, 0.98};