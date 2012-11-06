for(ix=0;ix<(snxy+1)*(snxy+1);ix++){
//for(iy=0;iy<=snxy;iy++){//y
int kk,iiy,iix,il;iix=(int)ix/(snxy+1);iiy=ix %(snxy+1);
for(kk=kmin;kk<=kmax;kk++){
doub t,maxy=fact*sftab[kk][1], xg=-maxy+2.*maxy/snxy*doub(iix),yg=-maxy+2.*maxy/snxy*doub(iiy),//*(1-1e-10)
b=sqrt(xg*xg+yg*yg),beta=atan2(yg,xg);//geon++;printf(" %d ",geon);//int currth = omp_get_thread_num();
#include "geoint.cpp"
ppy[currth].nu=1e9*sftab[kk][0];
#include "solvetrans.cpp"

(*ausin)[iix][iiy][kk][0]=II[0];
(*ausin)[iix][iiy][kk][1]=II[1]*cos(II[2])*sin(II[3]);
(*ausin)[iix][iiy][kk][2]=II[1]*sin(II[2])*sin(II[3]);
(*ausin)[iix][iiy][kk][3]=II[1]*cos(II[3]);
(*ausin)[iix][iiy][kk][4]=II[4];
};};//};

for(kk=kmin;kk<=kmax;kk++){
for(il=0;il<5;il++)in[kk][il]=0.;doub hei=2./snxy;doub maxy=fact*sftab[kk][1];

for(ix=0;ix<=snxy-3;ix+=2)for(iy=0;iy<=snxy-3;iy+=2)for(il=0;il<5;il++)
in[kk][il]+=((*ausin)[ix][iy][kk][il] + (*ausin)[ix][2 + iy][kk][il] + 4*(*ausin)[1 +ix][iy][kk][il] + 4*(*ausin)[1 + ix][2 + iy][kk][il] + (*ausin)[2 +ix][iy][kk][il] +
4*((*ausin)[ix][1 + iy][kk][il] + 4*(*ausin)[1 + ix][1 +iy][kk][il] + (*ausin)[2 + ix][1 + iy][kk][il]) + (*ausin)[2 + ix][2 +iy][kk][il])*hei*hei/9.;
for(ix=0;ix<=snxy-1;ix++)for(il=0;il<5;il++)in[kk][il]+=((*ausin)[ix][snxy-1][kk][il] + (*ausin)[ix][snxy][kk][il] + (*ausin)[1 + ix][snxy-1][kk][il] + (*ausin)[1 + ix][snxy][kk][il])*hei*hei/4.;
for(iy=0;iy<=snxy-2;iy++)for(il=0;il<5;il++)in[kk][il]+=((*ausin)[snxy-1][iy][kk][il] + (*ausin)[snxy-1][1 + iy][kk][il] + (*ausin)[snxy][iy][kk][il] + (*ausin)[snxy][1 + iy][kk][il])*hei*hei/4.;

for(il=0;il<5;il++)in[kk][il]*=maxy*maxy;
totin[kk]=66.4648*in[kk][0];LPo[kk]=100.*sqrt(in[kk][1]*in[kk][1]+in[kk][2]*in[kk][2])/in[kk][0];CP[kk]=100.*in[kk][3]/in[kk][0];
ang[kk]=fmod(180/PI*atan2(in[kk][2],in[kk][1])/2.+180.,180.);err[kk]=66.4648*in[kk][4];//err[kk]=100*fabs(in[kk][4]/in[kk][0]-1.);
printf("%d; f=%.1f; I=%.3fJy LP=%.2f%% EVPA=%.1fdeg CP=%.3f%% non-pol I=%.2fJy \n",fnum,sftab[kk][0], totin[kk],LPo[kk],ang[kk], CP[kk],err[kk]);
//Q-1;  U-2

};
xisq=pow((doub)(totin[4]-tofit[4][1])/dFnu[4],(doub)2.)+pow((doub)(totin[5]-tofit[5][1])/dFnu[5],(doub)2.)+pow((doub)(totin[6]-tofit[6][1])/dFnu[6],(doub)2.)+
pow((doub)(totin[7]-tofit[7][1])/dFnu[7],(doub)2.)+pow((doub)(totin[8]-tofit[8][1])/dFnu[8],(doub)2.)+pow((doub)(totin[9]-tofit[9][1])/dFnu[9],(doub)2.)+
pow((doub)(totin[10]-tofit[10][1])/dFnu[10],(doub)2.);dof=7.;
doub angshi=(-(ang[4]+ang[7]+ang[8])+(tofit[4][3]+tofit[7][3]+tofit[8][3]))/3.;xisq/=dof;

/*if(oo==0)printf("Fitting 88-657GHz spectrum: xi^2/dof=%.2f",xisq/dof);
if(wfit[0]){xisq+=pow((doub)(LPo[4]/tofit[4][2]-1)/dLP,(doub)2.)+pow((doub)(LPo[7]/tofit[7][2]-1)/dLP,(doub)2.)+pow((doub)(LPo[8]/tofit[8][2]-1)/dLP,(doub)2.);dof+=3.;};
if(wfit[1]){xisq+=pow((doub)(ang[4]+angshi-tofit[4][3])/dEVPA[0],(doub)2.)+pow((doub)(ang[7]+angshi-tofit[7][3])/dEVPA[1],(doub)2.)+pow((doub)(ang[8]+angshi-tofit[8][3])/dEVPA[2],(doub)2.),dof+=3.;};
if(wfit[2]){xisq+=pow((doub)(CP[7]-tofit[7][4])/dCP,(doub)2.)+pow((doub)(CP[8]-tofit[8][4])/dCP,(doub)2.);dof+=2.;};
*/
//if(oo==0){printf("; full xi^2/dof=%.2f; ang0=%.2f\n",xisq,angshi);xxisq=xisq;xang=angshi;};

for(il=4;il<=10;il++)resid[il-4][oo]=(totin[il]-tofit[il][1])/dFnu[il];//moving into residual vector
/*if(wfit[0]){resid[7][oo]=(LPo[4]/tofit[4][2]-1)/dLP;resid[8][oo]=(LPo[7]/tofit[7][2]-1)/dLP; resid[9][oo]=(LPo[8]/tofit[8][2]-1)/dLP;} else {for(il=7;il<=9;il++)resid[il][oo]=0.;};
if(wfit[1]){resid[10][oo]=(ang[4]+angshi-tofit[4][3])/dEVPA[0];resid[11][oo]=(ang[7]+angshi-tofit[7][3])/dEVPA[1]; resid[12][oo]=(ang[8]+angshi-tofit[8][3])/dEVPA[2];}else {for(il=10;il<=12;il++)resid[il][oo]=0.;};
if(wfit[2]){resid[13][oo]=(CP[7]-tofit[7][4])/dCP;resid[14][oo]=(CP[8]-tofit[8][4])/dCP;}else {for(il=13;il<=14;il++)resid[il][oo]=0.;}
*/
if(iswrite){
stringstream sstr;sstr <<(int)100*a<<"th"<<(int)floor(100*xth+1e-6)<<"fn"<<fnum<<"boo";//"boo" is a legacy string
string stra = sstr.str();
FILE * pFile; pFile=fopen ((dir+"poliresa"+stra+fif+".dat").c_str(),"a");
for(kk=kmin;kk<=kmax;kk++)
fprintf(pFile,"%d %.2f %.5f %.3f %.3f %.3f %.3f %.5f %.4f %.4f %.2f %.4e %.4f\n",
		fnum,sftab[kk][0],totin[kk],LPo[kk],ang[kk],CP[kk],err[kk],heat,rhonor,xxisq,TpTe,rate*year/Msun,th);
fclose(pFile);};