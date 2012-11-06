//int i;//double y4sq=y[4]*y[4], y5sq=y[5]*y[5], y45sq=y4sq + asq*y5sq, y45sqsq=y45sq*y45sq;
doub det=-3.*cosh(llog),t;
t=1.+ 3.* sinh(llog)/r0;//evalpoint(t);
//#include "evalpoint.cpp"
#include "evalpointxx.cpp"

doub cosqs=cos(yyy[2]),sinqs=sin(yyy[2]),frsq=fr*fr,sinths=sin(yyy[3]),cosths=cos(yyy[3]);
ff[0]=det*(jIc/frsq + fr*(-(aIc*yyy[0]) - aVc*cosths*yyy[1] - aQc*cos2k*cosqs*sinths*yyy[1] - aQc*sin2k*sinqs*sinths*yyy[1]));
ff[1]=det*((cosths*jVc + jQc*(cos2k*cosqs + sin2k*sinqs)*sinths)/frsq + fr*(-(aVc*cosths*yyy[0]) - aQc*(cos2k*cosqs + sin2k*sinqs)*sinths*yyy[0] - aIc*yyy[1]));
ff[2]=det*((jQc*(cosqs*sin2k - cos2k*sinqs))/(frsq*sinths*yyy[1]) + 
    (fr*(-((aQc*(cosqs*sin2k - cos2k*sinqs)*yyy[0])/sinths) + (rVc - (cos2k*cosqs*cosths*rQc)/sinths - (cosths*rQc*sin2k*sinqs)/sinths)*yyy[1]))/yyy[1]);
ff[3]=det*((cosths*jQc*(cos2k*cosqs + sin2k*sinqs) - jVc*sinths)/(frsq*yyy[1]) + 
    (fr*(-(aQc*cosths*(cos2k*cosqs + sin2k*sinqs)*yyy[0]) + aVc*sinths*yyy[0] + rQc*(cosqs*sin2k - cos2k*sinqs)*yyy[1]))/yyy[1]);
ff[4]= det*(-fr*aIc*yyy[4] + jIc/frsq);
/*
ff[0]=det*(jIc/frsq + fr*(-(aIc*yyy[0]) - aQc*cos2k*cosqs*yyy[1] - aQc*sin2k*sinqs*yyy[1] - aVc*yyy[3]));
ff[1]=det*((cos2k*cosqs*jQc + jQc*sin2k*sinqs)/frsq + fr*(-(aIc*yyy[1]) - aQc*cos2k*cosqs*yyy[0] + cosqs*rQc*sin2k*yyy[3] - sinqs*(aQc*sin2k*yyy[0] + cos2k*rQc*yyy[3])));
ff[2]=det*((cosqs*jQc*sin2k - cos2k*jQc*sinqs)/(yyy[1]*frsq) + 
   fr*(rVc*yyy[1] + aQc*cos2k*sinqs*yyy[0] - rQc*sin2k*sinqs*yyy[3] - cosqs*(aQc*sin2k*yyy[0] + cos2k*rQc*yyy[3]))/yyy[1]);

ff[3]=det*(jVc/frsq + fr*(-(aVc*yyy[0]) - cosqs*rQc*sin2k*yyy[1] + cos2k*rQc*sinqs*yyy[1] - aIc*yyy[3]));
ff[4]= det*(-fr*aIc*yyy[4] + jIc/frsq);

*/
//ff[0]=det*(jIc/fr/fr + fr*(-aIc*yyy[0] - aQc*cos2k*yyy[1] - aQc*sin2k*yyy[2] - aVc*yyy[3]));
//ff[1]= det*((cos2k*jQc)/fr/fr + fr*(-(aQc*cos2k*yyy[0]) - aIc*yyy[1] - rVc*yyy[2] + rQc*sin2k*yyy[3]));
//ff[2]= det*((jQc*sin2k)/fr/fr + fr*(-(aQc*sin2k*yyy[0]) + rVc*yyy[1] - aIc*yyy[2] - cos2k*rQc*yyy[3])); 
//ff[3]= det*(jVc/fr/fr + fr*(-(aVc*yyy[0]) - rQc*sin2k*yyy[1] + cos2k*rQc*yyy[2] - aIc*yyy[3]));
if(ff[0]!=ff[0])
{t=0.;}

return GSL_SUCCESS;
