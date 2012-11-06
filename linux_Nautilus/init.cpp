//using namespace std;
int m,k,j,tn,xrn,rn,n,sn,i,nr,nt,np;
bool fl=true;filebuf *pbuf;
doub rr, costh, temp,tempsq, cossq,sinsq,rsq,Del,rhosq;
doub iKS[4][4], iMKS[4][4], MKStoKS[4][4], rest[11], Bup[4],Bi[4],u[4],uKS[4],uloKS[4],ulo[4],BupKS[4],rho,sc[4];


a=atab[sp];asq=a*a;ncut=ncuttab[sp];
fdir=adir+astr[sp]+fieldstr;
if(fdiff>fdiffmax){printf("Can't load that many snapshots... \n");return 0;};
mintim=1-0.0001*fdiff;maxtim=1+0.0001*fdiff;

for(i=-fdiff;i<=fdiff;i++){
	stringstream sstr;
	sstr<<setfill('0')<<setw(4)<<fnum+i;
	ifstream fline((fdir+"fieldline"+sstr.str()+".bin").c_str(),ios::in|ios::binary);
	pbuf=fline.rdbuf();
	int fsize=pbuf->pubseekoff (0,ios::end),tosize=11*phlen*thlen*rlen*sizeof(float);
	pbuf->pubseekpos(fsize-tosize);
	fline.read(reinterpret_cast<char *>((*uu)[fdiff+i]), tosize);
if(fline.good())printf("fieldline %d loaded\n",fnum+i); else {printf("fieldline %d cannot be loaded \n",fnum+i);exit(1);}
	fline.close();}
//rho, u, -u^t, -T^r_t/(rho u^r), u^t, v^r, v^theta, v^phi, B^r, B^theta, B^phi

if(!inited){
//r, theta, rho (density), U (internal energy) , u^0, u^1, u^2, u^3, u_0, u_1, u_2, u_3, b_1, b_2, b_3, (b_1)^2, (b_2)^2, (b_3)^2
	//reading or creating density & temperature profiles
	ifstream faire ((dir+astr[sp]+xstr+"Tsmap"+astr[sp]+".dat").c_str(), ios::in);
	for(k=0;k<rlen;k++)for(i=0;i<dd;i++)faire>>xx[k][i];
if(faire.good())printf("Tsmap loaded\n"); else {printf("Tsmap cannot be loaded \n");exit(2);}
	faire.close();

//	ifstream gre((dir+astr[sp]+xstr+"usgdump2d").c_str(), ios::in|ios::binary);
//	pbuf=gre.rdbuf(); pbuf->pubseekpos(usgoff);int tosize=usgsize*rlen*thlen*sizeof(double);
//	gre.read(reinterpret_cast<char *>(*usgread), tosize);
//	gre.close();

	ifstream dxp((dir+astr[sp]+xstr+"dxdxp.dat").c_str(), ios::in|ios::binary);
	pbuf=dxp.rdbuf();
	dxp.read(reinterpret_cast<char *>(coord), ndd*thlen*2*sizeof(float));
	dxp.read(reinterpret_cast<char *>(dxdxp), ndd*thlen*4*4*sizeof(float));
	dxp.close();

	for(k=0;k<ndd;k++)for(i=0;i<thlen;i++){
	//theta[k][i]=-cos((*usgread)[i][k][8]);//reading any grid
	theta[k][i]=coord[k][i][1];//reading any grid
	};

	for(k=0;k<ncut;k++){
	//rtab[k]=(*usgread)[0][k][7];
	rtab[k]=exp(coord[k][0][0]);
	Tstab[k]=xx[k][1];};maxT=Tstab[0];//reading radial grid

	ifstream fI ((dir+"lookupjIy.dat").c_str(), ios::in);
	ifstream fQ ((dir+"lookupjQy.dat").c_str(), ios::in);
	ifstream fV ((dir+"lookupjVy.dat").c_str(), ios::in);
	ifstream frQ ((dir+"lookuprQa.dat").c_str(), ios::in);
	ifstream frV ((dir+"lookuprVa.dat").c_str(), ios::in);
	for(k=0;k<=Tlen;k++) for(j=0;j<=nWlen;j++) fI>>jI[k][j];
	for(k=0;k<=Tlen;k++) for(j=0;j<=nWlen;j++) fQ>>jQ[k][j];
	for(k=0;k<=Tlen;k++) for(j=0;j<=nWlen;j++) fV>>jV[k][j];
	for(k=0;k<=2*Tlen;k++) frQ>>rQ[k];
	for(k=0;k<=2*Tlen;k++) frV>>rV[k];
	//for(k=0;k<=Tlen;k++) for(j=0;j<=nWlen;j++) frQ>>rQ[k][j];
	//for(k=0;k<=Tlen;k++) for(j=0;j<=nWlen;j++) frV>>rV[k][j];
	fI.close();fQ.close();fV.close();frQ.close();frV.close();

	rmin=xx[0][0]; rmax=xx[rlen-1][0];
	//lrmin=log(rmin-off);lrmax=log(rmax-off);
	inited=true;};

rcut=xx[ncut-1][0];//switching to realistic boundary at 3.4*10^5M
rhopo=-log(rhoout*dense/rhonor/xx[ncut-1][2])/log(rrmax/rcut);
Upo=-log(Tout/xx[ncut-1][1])/log(rrmax/rcut);
// for Bpo=2.0 => extension of radial profile is B~1/r, what makes plasma beta approximately constant
//----------------tests------------------------
	if(echeck1)rhopo+=0.2;
	if(echeck2)Upo-=0.1;
	if(echeck3)Bpo+=0.2;
//--------------end---tests--------------------

rhocon=rhonor*xx[ncut-1][2];
Ucon=xx[ncut-1][1];
Bnor=sqrt(4.*PI*rhonor*mp)*cc; //can set "-" sign

for(k=ncut;k<ndd;k++){rtab[k]=rcut*pow((doub)rrmax/rcut,(k-ncut+(doub)1.)/(doub)(ndd-ncut));//extending radial grid
Tstab[k]=Ucon*pow((doub)rtab[k]/rcut,(doub)-Upo);};
//for(k=0;k<ndd;k++)for(i=0;i<thlen;i++){coord[k][i][0]=log(rtab[k]);coord[k][i][1]=theta[k][i];};//filling in the entire coordinate grid
int nx=60;doub rx=exp(coord[nx-1][0][0]),u0;
rate=0.;
for(k=0;k<thlen-1;k++)for(i=0;i<phlen;i++)
	for(j=0;j<=2*fdiff;j++) {u0=(*uu)[j][i][k][nx-1][4];usp[j][i][k][0]=u0; usp[j][i][k][1]=u0*(*uu)[j][i][k][nx-1][5];usp[j][i][k][2]=u0*(*uu)[j][i][k][nx-1][6];usp[j][i][k][3]=u0*(*uu)[j][i][k][nx-1][7];
	uspKS[j][i][k]=0.;for(m=0;m<4;m++)uspKS[j][i][k]+=dxdxp[nx-1][k][1][m]*usp[j][i][k][m];//radial velocity in KS metric
rate+=2./phlen*PI*(rx*rx+ asq*theta[nx-1][k]*theta[nx-1][k+1])*(*uu)[j][i][k][nx-1][0]*uspKS[j][i][k]*(theta[nx-1][k]-theta[nx-1][k+1]);};
rate*=rhonor*rgrav*rgrav*cc*mp/(2*fdiff+1);

//int nx=60,m;doub rx=exp(coord[nx-1][0][0]);rate=0.;
//for(k=0;k<thlen-1;k++)for(i=0;i<phlen;i++)
//	for(j=0;j<=2*fdiff;j++) rate+=2./phlen*PI*(rx-off)*(rx*rx+ asq*theta[nx-1][k]*theta[nx-1][k+1])*(*uu)[j][i][k][nx-1][0]*(*uu)[j][i][k][nx-1][4]*(*uu)[j][i][k][nx-1][5]*(theta[nx-1][k]-theta[nx-1][k+1]);
//rate*=rhonor*rgrav*rgrav*cc*mp/(2*fdiff+1);

doub acc=1e-5, IT[2], hz,rz, step;//int k; //abs accuracy 1e-6 and method 
gsl_odeiv_step *sz; gsl_odeiv_control *cz; gsl_odeiv_evolve *ez;Tz=gsl_odeiv_step_rk2;
gsl_odeiv_system sysT = {solte, NULL, 2,NULL};
ez = gsl_odeiv_evolve_alloc(2);
cz = gsl_odeiv_control_standard_new(0.0, acc, 1.0, 0.0);
sz = gsl_odeiv_step_alloc (Tz, 2);

IT[0]=Tstab[ndd-1];IT[1]=IT[0];stNt=0;rz=rrmax;
step=-0.0025*rz;hz=step;
ts[0]=IT[0];tp[0]=IT[0];te[0]=IT[0];

while (rz > 1.001*rmin)
{	stNt++; if(-hz>0.012*rz)hz=-0.008*rz;
	int status = gsl_odeiv_evolve_apply (ez, cz, sz, &sysT, &rz, 1.001*rmin, &hz, IT);
te[stNt]=IT[0];
tp[stNt]=IT[1];
if((rz<6) && (fl)){fl=false;TpTe=tp[stNt]/te[stNt];Te6=te[stNt];
printf("fn=%d stN=%d r=%.2fM Tp/Te=%.2f Te=%.3e rate=%.3e he=%.3f rho=%.3e\n", fnum, stNt,rz, TpTe, Te6, rate*year/Msun,heat,rhonor);};
;};
maxT=ts[stNt];minT=ts[0];

if(stNt>10000)printf("Exceeded buffer! Recalculate!\n");
gsl_odeiv_evolve_free (ez);gsl_odeiv_control_free (cz);gsl_odeiv_step_free (sz);
//compute quantities in the co-moving frame at the boundary
xrn=ncut-1;sn=fdiff;
rr=rcut;//=xx[ncut-1][0];

for(nt=0;nt<thlen;nt++)//for(nr=0;nr<rlen;nr++)
	{costh=theta[xrn][nt];
	cossq=costh*costh;sinsq=1-cossq;rsq=rr*rr;rhosq=rsq+asq*cossq;Del=rsq-2.*rr+asq;
	doub temp=2.*rr/rhosq;//-,+,+,+ signature covariant KS matrix
	iKS[0][0]=temp-1.;iKS[0][1]=temp;iKS[0][2]=0.;iKS[0][3]=-a*temp*sinsq;
	iKS[1][0]=iKS[0][1]; iKS[1][1]=1.+temp;iKS[1][2]=0.;iKS[1][3]=-a*(1.+temp)*sinsq;
	iKS[2][0]=iKS[0][2];iKS[2][1]=iKS[1][2];iKS[2][2]=rhosq;iKS[2][3]=0.;
	iKS[3][0]=iKS[0][3];iKS[3][1]=iKS[1][3];iKS[3][2]=iKS[2][3];iKS[3][3]=sinsq*(rhosq+asq*(1+temp)*sinsq);

	for(i=0;i<4;i++)for(k=0;k<4;k++)MKStoKS[i][k]=dxdxp[xrn][nt][i][k];
	for(i=0;i<4;i++)for(k=0;k<4;k++)iMKS[i][k]=0.;
	for(i=0;i<4;i++)for(k=0;k<4;k++)for(j=0;j<4;j++)for(m=0;m<4;m++)iMKS[i][k]+=MKStoKS[m][i]*iKS[m][j]*MKStoKS[j][k];
	
	for(np=0;np<phlen;np++){
		for(m=0;m<11;m++)rest[m]=(*uu)[sn][np][nt][xrn][m];
		rho=rest[0]*rhonor; rest[1]*=mp*cc*cc/3/kb/rest[0];
		u[0]=rest[4];u[1]=u[0]*rest[5];u[2]=u[0]*rest[6];u[3]=u[0]*rest[7];
		Bi[1]=Bnor*rest[8];Bi[2]=Bnor*rest[9];Bi[3]=Bnor*rest[10];
uext[np][nt][3]=rho;
uext[np][nt][4]=rest[1];

		for(i=1;i<4;i++){uKS[i]=0.;for(k=0;k<4;k++)uKS[i]+=MKStoKS[i][k]*u[k];}
		doub expr=fabs((uKS[1]*iKS[0][1] + uKS[3]*iKS[0][3])* (uKS[1]*iKS[0][1] + uKS[3]*iKS[0][3]) - iKS[0][0]*(1. + uKS[1]*uKS[1]*iKS[1][1] + 2.*uKS[1]*uKS[3]*iKS[1][3] + uKS[2]*uKS[2]*iKS[2][2] + uKS[3]*uKS[3]*iKS[3][3]));
		uKS[0]=-((uKS[1]*iKS[0][1] + uKS[3]*iKS[0][3] + sqrt(expr))/iKS[0][0]);
		for (m=0;m<4;m++){uloKS[m]=0.;for(j=0;j<4;j++)uloKS[m]+=iKS[m][j]*uKS[j];};//uloKS
		u[0]=uKS[0];//works in simple MKS, where only spatial transformation is performed
		for(m=0;m<4;m++){ulo[m]=0.;for(j=0;j<4;j++)ulo[m]+=iMKS[m][j]*u[j];};//ulo

		doub yyKS[4][4]={{uKS[0], uKS[1], uKS[2], uKS[3]},{uKS[0]*uloKS[1], -(uKS[0]*uloKS[0] + uKS[3]*uloKS[3]), 0, uKS[3]*uloKS[1]},
		{uloKS[2]*uKS[0], uloKS[2]*uKS[1],uloKS[2]*uKS[2]+1, uloKS[2]*uKS[3]},{-uloKS[3]/uloKS[0], 0, 0, 1}}, yyloKS[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},
		xprx[4]={uKS[0]*uloKS[1], -(uKS[0]*uloKS[0] + uKS[3]*uloKS[3]), 0, uKS[3]*uloKS[1]},xloprx[4]={0,0,0,0};
		for(j=0;j<4;j++)for(n=0;n<4;n++)xloprx[j]+=iKS[j][n]*xprx[n];
		for(m=0;m<4;m++)for(j=0;j<4;j++)for(n=0;n<4;n++)yyloKS[m][n]+=iKS[j][n]*yyKS[m][j];//already transposed matrix

		for (m=0;m<4;m++){sc[m]=0.;for(j=0;j<4;j++)sc[m]+=yyloKS[m][j]*yyKS[m][j];
		sc[m]=sqrt(fabs(sc[m]));for(j=0;j<4;j++){yyloKS[m][j]/=sc[m];yyKS[m][j]/=sc[m];};};
		for(j=0;j<4;j++)yyloKS[0][j]=-yyloKS[0][j];
		for(m=0;m<4;m++){yyKS[1][m]=xprx[m];for(j=0;j<4;j++)yyKS[1][m]-=yyKS[3][m]*xloprx[j]*yyKS[3][j];};
		//updating x1 vector;

		for(j=0;j<4;j++){yyloKS[1][j]=0.;for(n=0;n<4;n++)yyloKS[1][j]+=iKS[n][j]*yyKS[1][n];};
		sc[1]=0.;for(j=0;j<4;j++)sc[1]+=yyloKS[1][j]*yyKS[1][j];sc[1]=sqrt(fabs(sc[1]));
		for(j=0;j<4;j++){yyloKS[1][j]/=sc[1];yyKS[1][j]/=sc[1];}

		doub blo=(ulo[1]*Bi[1]+ulo[2]*Bi[2]+ulo[3]*Bi[3]);
		Bup[0]=blo;for(m=1;m<4;m++)Bup[m]=(Bi[m]+blo*u[m])/u[0];
		for(i=0;i<4;i++){BupKS[i]=0.;for(k=0;k<4;k++)BupKS[i]+=MKStoKS[i][k]*Bup[k];}
		doub test[4]={0.,0.,0.,0.};for(m=0;m<4;m++){for(j=0;j<4;j++)test[m]+=yyloKS[m][j]*uKS[j];};
		doub B[4]={0.,0.,0.,0.};for(m=0;m<4;m++){for(j=0;j<4;j++)B[m]+=yyloKS[m][j]*BupKS[j];};

uext[np][nt][0]=B[1];
uext[np][nt][1]=B[2];
uext[np][nt][2]=B[3];
		};
};
printf("\n");return(0);
//rest[m]=(*uu)[sn][k][tn][xrn][m]