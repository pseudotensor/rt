//int i;//double y4sq=y[4]*y[4], y5sq=y[5]*y[5], y45sq=y4sq + asq*y5sq, y45sqsq=y45sq*y45sq;
doub det=-3.*cosh(llog),t;
t=1.+ 3.* sinh(llog)/r0;//evalpoint(t);
//#include "evalpoint.cpp"
#include "evalpointxx.cpp"

ff[0]=det*(jIc/fr/fr + fr*(-aIc*yyy[0] - aQc*cos2k*yyy[1] - aQc*sin2k*yyy[2] - aVc*yyy[3]));
ff[1]= det*((cos2k*jQc)/fr/fr + fr*(-(aQc*cos2k*yyy[0]) - aIc*yyy[1] - rVc*yyy[2] + rQc*sin2k*yyy[3]));
ff[2]= det*((jQc*sin2k)/fr/fr + fr*(-(aQc*sin2k*yyy[0]) + rVc*yyy[1] - aIc*yyy[2] - cos2k*rQc*yyy[3])); 
ff[3]= det*(jVc/fr/fr + fr*(-(aVc*yyy[0]) - rQc*sin2k*yyy[1] + cos2k*rQc*yyy[2] - aIc*yyy[3]));
ff[4]= det*(-fr*aIc*yyy[4] + jIc/fr/fr);
if(ff[0]!=ff[0])
{t=0.;}

return GSL_SUCCESS;