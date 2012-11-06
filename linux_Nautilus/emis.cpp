//double emis (double nW, double tet)
doub lnW,lTe,lTefrac,lTeman,lnWfrac,lnWman;
lnW=log(nW);lTe=log(kb*tet/me/cc/cc);

if((lTe<lTmin)||(lTe>lTmax)||(lnW<lnWmin)||(lnW>lnWmax)){xVc=0.;xQc=0.;xIc=0.;} else {
lTefrac=(lTe-lTmin)/(lTmax-lTmin)*Tlen;j=floor(lTefrac);lTeman=lTefrac-j;//x2
lnWfrac=(lnW-lnWmin)/(lnWmax-lnWmin)*nWlen;k=floor(lnWfrac);lnWman=lnWfrac-k;//lr
temp = lnWman*lTeman;
xIc=exp((-lnWman - lTeman + 1 + temp)*jI[j][k] + (lnWman - temp)*jI[j][1 + k] + (lTeman - temp)*jI[1 + j][k] + temp*jI[1 + j][1 + k]);
xQc=exp((-lnWman - lTeman + 1 + temp)*jQ[j][k] + (lnWman - temp)*jQ[j][1 + k] + (lTeman - temp)*jQ[1 + j][k] + temp*jQ[1 + j][1 + k]);
xVc=exp((-lnWman - lTeman + 1 + temp)*jV[j][k] + (lnWman - temp)*jV[j][1 + k] + (lTeman - temp)*jV[1 + j][k] + temp*jV[1 + j][1 + k]);};

if((lTe<lTminr)||(lTe>lTmax)){xrQc=0.;xrVc=0.;} else {
lTefrac=(lTe-lTminr)/(lTmax-lTminr)*2*Tlen;j=floor(lTefrac);lTeman=lTefrac-j;//x2
xrQc=exp((-lTeman + 1)*rQ[j] + (lTeman)*rQ[1 + j]);
xrVc=exp((-lTeman + 1)*rV[j] + (lTeman)*rV[1 + j]);};

if(xIc!=xIc){
m=1;}
//to be improved - 