{int ind=co,fmin=6734,fmax=7034;sep=(fmax-fmin)/(ind-1);FILE * yfline;filebuf *pbuf;
rhonor=1000000.; heat=0.5;  th=1.0;fdiff=0;
	int m,k,j,n,z, nx1,nx2,nx3,rk,thk,phk;
doub h1, dx1,x1min,tt,off,x1ar[rlen],rad[rlen],rho[rlen],ts[rlen];

	fdir=adir+astr[sp]+fieldstr;a=atab[sp];asq=a*a;
	stringstream sstr;
	sstr<<setfill('0')<<setw(4)<<fmin;
	ifstream xfline((fdir+"fieldline"+sstr.str()+".bin").c_str(),ios::in);
	xfline>>tt>>nx1>>nx2>>nx3>>x1min>>dx1>>dx1>>dx1;
	xfline>>off>>off>>off>>off>>off>>off;
	xfline>>h1>>h1>>h1;xfline.close();
if((rlen!=nx1)||(thlen!=nx2)||(phlen!=nx3)){
printf("Errors in dimensions");return 0;};


//a=atab[sp];asq=a*a;ncut=ncuttab[sp];
//fdir=adir+astr[sp]+fieldstr;
if(astr[sp].length()<4){
	ifstream gre((dir+astr[sp]+xstr+"usgdump2d").c_str(), ios::in|ios::binary);
	pbuf=gre.rdbuf(); pbuf->pubseekpos(usgoff);int tosize=usgsize*rlen*thlen*sizeof(double);
	gre.read(reinterpret_cast<char *>(*usgread), tosize);
	gre.close();
for(i=0;i<rlen;i++){rad[i]=(*usgread)[0][i][7];}//creating radial grid
	for(k=0;k<rlen;k++)for(i=0;i<thlen;i++)theta[k][i]=-cos((*usgread)[i][k][8]);//reading theta grid
} else
{
	ifstream dxp((dir+astr[sp]+xstr+"dxdxp.dat").c_str(), ios::in|ios::binary);
	pbuf=dxp.rdbuf();
	dxp.read(reinterpret_cast<char *>(coord), ndd*thlen*2*sizeof(float));
	dxp.read(reinterpret_cast<char *>(dxdxp), ndd*thlen*4*4*sizeof(float));
	dxp.close();
	for(k=0;k<rlen;k++)for(i=0;i<thlen;i++){theta[k][i]=coord[k][i][1];}//reading any grid
	for(k=0;k<rlen;k++){rad[k]=exp(coord[k][0][0]);}
}
for(rk=0;rk<rlen;rk++){rho[rk]=0.;ts[rk]=0.;}

	for(fnum=fmin;fnum<=fmax;fnum+=sep){
	stringstream sstr;
	sstr<<setfill('0')<<setw(4)<<fnum;
	ifstream fline((fdir+"fieldline"+sstr.str()+".bin").c_str(),ios::in|ios::binary);
	pbuf=fline.rdbuf();
	int fsize=pbuf->pubseekoff (0,ios::end),tosize=11*phlen*thlen*rlen*sizeof(float);
	pbuf->pubseekpos(fsize-tosize);
	fline.read(reinterpret_cast<char *>(uu[0]), tosize);
	fline.close();
for(rk=0;rk<rlen;rk++)for(thk=thlen/2-3;thk<=thlen/2+3;thk++)for(phk=0;phk<phlen;phk++)
{rho[rk]+=(*uu[0])[phk][thk][rk][0]/7/phlen/ind;
ts[rk]+=(*uu[0])[phk][thk][rk][1]/7/phlen/ind*mp*cc*cc/3/kb/(*uu[0])[phk][thk][rk][0];};

//int nx=37;doub rx=rad[nx-1]; //for(m=20;m<=200;m+=10){nx=m;rx=xx[nx-1][0];
	//rate=0.;
//	for(k=0;k<thlen-1;k++)for(i=0;i<phlen;i++)
//	rate+=2./phlen*PI*(rx-off)*(rx*rx+ asq*theta[0][k]*theta[0][k+1])*(*uu[0])[i][k][nx-1][0]*(*uu[0])[i][k][nx-1][4]*(*uu[0])[i][k][nx-1][5]*(theta[0][k]-theta[0][k+1]);
//	rate*=rhonor*rgrav*rgrav*cc*mp;
//


int nx=37,m;doub rx=rad[nx-1],u0;//can't compute rate without knowing dxdxp, for which another calculations is required
//rate=0.;
//for(k=0;k<thlen-1;k++)for(i=0;i<phlen;i++){u0=(*uu[0])[i][k][nx-1][4];usp[0][i][k][0]=u0; usp[0][i][k][1]=u0*(*uu[0])[i][k][nx-1][5];usp[0][i][k][2]=u0*(*uu[0])[i][k][nx-1][6];usp[0][i][k][3]=u0*(*uu[0])[i][k][nx-1][7];
//	uspKS[0][i][k]=0.;for(m=0;m<4;m++)uspKS[0][i][k]+=dxdxp[nx-1][k][1][m]*usp[0][i][k][m];//radial velocity in KS metric
//rate+=2./phlen*PI*(rx*rx+ asq*theta[nx-1][k]*theta[nx-1][k+1])*(*uu[0])[i][k][nx-1][0]*uspKS[0][i][k]*(theta[nx-1][k]-theta[nx-1][k+1]);};
//rate*=rhonor*rgrav*rgrav*cc*mp;
//printf("fnum=%d  rate=%.3e\n", fnum, rate*year/Msun);
doub rhoinst=0.;
for(thk=thlen/2-3;thk<=thlen/2+3;thk++)for(phk=0;phk<phlen;phk++)rhoinst+=(*uu[0])[phk][thk][nx-1][0]/7/phlen;
printf("fn=%d: average density at r=%f is rho=%f \n",fnum,rx,rhoinst);
}
	

bool sorted=false;doub temp;
while(!sorted){sorted=true;
for(ix=1;ix<rlen;ix++)
{if(ts[ix-1]<ts[ix]){temp=ts[ix-1];ts[ix-1]=ts[ix];ts[ix]=temp;sorted=false;}
if(rho[ix-1]<rho[ix]){temp=rho[ix-1];rho[ix-1]=rho[ix];rho[ix]=temp;sorted=false;};}
;}

yfline=fopen((adir+astr[sp]+xstr+"Tsmapx"+astr[sp]+".dat").c_str(),"w");
for(rk=0;rk<rlen;rk++) 
	fprintf(yfline,"%.8e %.8e %.8e \n",rad[rk],ts[rk],rho[rk]);
fclose(yfline);

}