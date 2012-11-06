//int i;//double y4sq=y[4]*y[4], y5sq=y[5]*y[5], y45sq=y4sq + asq*y5sq, y45sqsq=y45sq*y45sq;
doub det=-3.*cosh(llog),t;
t=1.+ 3.* sinh(llog)/r0;
#include "evalpointzero.cpp"
//#include "evalpointnew.cpp"
//#include "evalpoint2012.cpp"

doub cosqs=cos(yyy[2]),sinqs=sin(yyy[2]),frsq=fr*fr, frcu=frsq*fr,sinths=sin(yyy[3]),cosths=cos(yyy[3]);

ff[0]=(det*(jIc - frcu*(aIc*yyy[0] + (aVc*cosths +aQc*cos2k*cosqs*sinths + aQc*sin2k*sinqs*sinths)*yyy[1])))/frsq;
ff[2]=(det*(cosqs*sin2k - cos2k*sinqs)*(jQc - aQc*frcu*yyy[0]) - 
det*frcu*(cos2k*cosqs*cosths*rQc + cosths*rQc*sin2k*sinqs -rVc*sinths)*yyy[1])/(frsq*sinths*yyy[1]);
ff[1]= (det*((cos2k*cosqs + sin2k*sinqs)*sinths*(jQc - aQc*frcu*yyy[0]) + cosths*(jVc -aVc*frcu*yyy[0]) - aIc*frcu*yyy[1]))/frsq;
ff[3]=(det*(-(jVc*sinths) + aVc*frcu*sinths*yyy[0] + cos2k*cosqs*cosths*(jQc - aQc*frcu*yyy[0]) + cosths*sin2k*sinqs*(jQc
- aQc*frcu*yyy[0]) + cosqs*frcu*rQc*sin2k*yyy[1] -cos2k*frcu*rQc*sinqs*yyy[1]))/(frsq*yyy[1]);
ff[4]= det*(-fr*xaIc*yyy[4] + xjIc/frsq);
//if(rr<10)ff[4]=-0.0001*(fr-1); else ff[4]=0.;
//if(rr<12)ff[4]=-0.001*u[3]; else ff[4]=0.;

if((fabs(ff[0])>1000.*fabs(ff[0]))&&(fabs(ff[0])>1e-6))
{i=1;}

//ff[0]=det*(jIc/frsq + fr*(-(aIc*yyy[0]) - aVc*cosths*yyy[1] - aQc*cos2k*cosqs*sinths*yyy[1] - aQc*sin2k*sinqs*sinths*yyy[1]));
//ff[1]=det*((cosths*jVc + jQc*(cos2k*cosqs + sin2k*sinqs)*sinths)/frsq + fr*(-(aVc*cosths*yyy[0]) - aQc*(cos2k*cosqs + sin2k*sinqs)*sinths*yyy[0] - aIc*yyy[1]));
//ff[2]=det*((jQc*(cosqs*sin2k - cos2k*sinqs))/(frsq*sinths*yyy[1]) + 
//    (fr*(-((aQc*(cosqs*sin2k - cos2k*sinqs)*yyy[0])/sinths) + (rVc - (cos2k*cosqs*cosths*rQc)/sinths - (cosths*rQc*sin2k*sinqs)/sinths)*yyy[1]))/yyy[1]);
//ff[3]=det*((cosths*jQc*(cos2k*cosqs + sin2k*sinqs) - jVc*sinths)/(frsq*yyy[1]) + 
//    (fr*(-(aQc*cosths*(cos2k*cosqs + sin2k*sinqs)*yyy[0]) + aVc*sinths*yyy[0] + rQc*(cosqs*sin2k - cos2k*sinqs)*yyy[1]))/yyy[1]);
//ff[4]= det*(-fr*aIc*yyy[4] + jIc/frsq);

//ff[0]=det*(jIc/fr/fr + fr*(-aIc*yyy[0] - aQc*cos2k*yyy[1] - aQc*sin2k*yyy[2] - aVc*yyy[3]));
//ff[1]= det*((cos2k*jQc)/fr/fr + fr*(-(aQc*cos2k*yyy[0]) - aIc*yyy[1] - rVc*yyy[2] + rQc*sin2k*yyy[3]));
//ff[2]= det*((jQc*sin2k)/fr/fr + fr*(-(aQc*sin2k*yyy[0]) + rVc*yyy[1] - aIc*yyy[2] - cos2k*rQc*yyy[3])); 
//ff[3]= det*(jVc/fr/fr + fr*(-(aVc*yyy[0]) - rQc*sin2k*yyy[1] + cos2k*rQc*yyy[2] - aIc*yyy[3]));
if(ff[0]!=ff[0])
{t=0.;printf("Error in radiative transfer\n");}
//if(yyy[0]<0){printf("haha %f, ix=%d, iy=%d",llog,0,0);}

return GSL_SUCCESS;
