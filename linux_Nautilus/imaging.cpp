for(ix=0;ix<(nxy+1)*(nxy+1);ix++){
//for(iy=0;iy<=nxy;iy++){//y
int kk,iiy,iix,il;iix=(int)ix/(nxy+1);iiy=ix %(nxy+1);
for(kk=kmin;kk<=kmax;kk++){
doub t,maxy=fact*sftab[kk][1], xg=-maxy+2.*maxy/nxy*doub(iix),yg=-maxy+2.*maxy/nxy*doub(iiy),//*(1-1e-10)
b=sqrt(xg*xg+yg*yg),beta=atan2(yg,xg);
#include "geoint.cpp"
ppy[currth].nu=1e9*sftab[kk][0];
#include "solvetrans.cpp"

/*(*intab)[ix][iy][0]=II[0];
(*intab)[ix][iy][1]=II[1]*cos(II[2])*sin(II[3]);
(*intab)[ix][iy][2]=II[1]*sin(II[2])*sin(II[3]);
(*intab)[ix][iy][3]=II[1]*cos(II[3]);
(*intab)[ix][iy][4]=II[4];*/

(*ausin)[iix][iiy][kk][0]=II[0];
(*ausin)[iix][iiy][kk][1]=II[1]*cos(II[2])*sin(II[3]);
(*ausin)[iix][iiy][kk][2]=II[1]*sin(II[2])*sin(II[3]);
(*ausin)[iix][iiy][kk][3]=II[1]*cos(II[3]);
(*ausin)[iix][iiy][kk][4]=ppy[currth].llmax;
//(*ausin)[ix][iy][kk][4]=II[4];

//if((35<iy)&&(iy<65)&&(85<ix)&&(ix<115))(*ausin)[ix][iy][kk][0]/=2.;
//if((125<iy)&&(iy<155)&&(85<ix)&&(ix<115))(*ausin)[ix][iy][kk][0]/=10.;
};};//};

for(kk=kmin;kk<=kmax;kk++){
for(ix=0;ix<=nxy;ix++)for(iy=0;iy<=nxy;iy++)for(il=0;il<=4;il++)(*intab)[ix][iy][il]=(*ausin)[ix][iy][kk][il];
for(il=0;il<5;il++)in[kk][il]=0.;doub hei=2./nxy;doub maxy=fact*sftab[kk][1];

for(ix=0;ix<=nxy-3;ix+=2)for(iy=0;iy<=nxy-3;iy+=2)for(il=0;il<=4;il++)
in[kk][il]+=((*intab)[ix][iy][il] + (*intab)[ix][2 + iy][il] + 4*(*intab)[1 +ix][iy][il] + 4*(*intab)[1 + ix][2 + iy][il] + (*intab)[2 +ix][iy][il] +
4*((*intab)[ix][1 + iy][il] + 4*(*intab)[1 + ix][1 +iy][il] + (*intab)[2 + ix][1 + iy][il]) + (*intab)[2 + ix][2 +iy][il])*hei*hei/9.;
for(ix=0;ix<=nxy-1;ix++)for(il=0;il<5;il++)in[kk][il]+=((*intab)[ix][nxy-1][il] + (*intab)[ix][nxy][il] + (*intab)[1 + ix][nxy-1][il] + (*intab)[1 + ix][nxy][il])*hei*hei/4.;
for(iy=0;iy<=nxy-2;iy++)for(il=0;il<5;il++)in[kk][il]+=((*intab)[nxy-1][iy][il] + (*intab)[nxy-1][1 + iy][il] + (*intab)[nxy][iy][il] + (*intab)[nxy][1 + iy][il])*hei*hei/4.;

for(il=0;il<5;il++)in[kk][il]*=maxy*maxy;
totin[kk]=66.4648*in[kk][0];LPo[kk]=100.*sqrt(in[kk][1]*in[kk][1]+in[kk][2]*in[kk][2])/in[kk][0];CP[kk]=100.*in[kk][3]/in[kk][0];
ang[kk]=fmod(180/PI*atan2(in[kk][2],in[kk][1])/2.+180.,180);err[kk]=66.4648*in[kk][4];//err[kk]=100*fabs(in[kk][4]/in[kk][0]-1.);
printf("%d; f=%.1f; I=%.3fJy LP=%.2f%% EVPA=%.1fdeg CP=%.3f%% non-pol I=%.2fJy \n",fnum,sftab[kk][0], totin[kk],LPo[kk],ang[kk], CP[kk],err[kk]);
//Q-1;  U-2


/*for(il=0;il<5;il++)in[il]*=maxy*maxy;
printf("I=%.3e, Q=%.3e, U=%.3e, V=%.3e, Ix=%.3e\n",in[0],in[1],in[2],in[3],in[4]);
printf("I=%.3fJy  LP=%.3f%%  CP=%.3f%%  non-pol error=%.3f%% \n",66.4648*in[0],100.*sqrt(in[1]*in[1]+in[2]*in[2])/in[0],100.*fabs(in[3])/in[0],100*fabs(in[4]/in[0]-1.));
*/

(*params)[0]=a;(*params)[1]=th;(*params)[2]=double(nxy);(*params)[3]=double(sftab[kk][0]);(*params)[4]=double(sftab[kk][1]);
(*params)[5]=double(heat);(*params)[6]=double(rhonor);
(*params)[7]=double(totin[kk]),(*params)[8]=double(LPo[kk]),(*params)[9]=double(ang[kk]),
(*params)[10]=double(CP[kk]);
(*params)[11]=double(err[kk]);(*params)[12]=double(TpTe);(*params)[13]=rate*year/Msun;

stringstream sstr;sstr<<(int)100*a<<"th"<<(int)floor(100*tth+1e-6)<<"f"<<floor(sftab[kk][0])<<"fn"<<fnum<<"_"<<nxy;string stra = sstr.str();

ofstream faire ((dir+"shotimag"+stra+".dat").c_str(), ios::out|ios::binary);
faire.write(reinterpret_cast<char *>(params), 20*sizeof(double));
faire.write(reinterpret_cast<char *>(intab), 5*(nxy+1)*(nxy+1)*sizeof(doub));
faire.close();
}
//WinExec("\"C:\\Program Files\\Wolfram Research\\Mathematica\\7.0\\math\" -noprompt -initfile \"G:\\mathematica\\3DMHD\\pol_imaging_script\"",SW_SHOWDEFAULT);