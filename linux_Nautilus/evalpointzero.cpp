bool ju=true,sing;
doub rr, costh, lrman,thman,ph,phfrac,phman, timfrac, timman, temp,tempsq, cossq,sinsq,rsq,Del,rhosq;
doub iKS[4][4], iMKS[4][4], MKStoKS[4][4], e1BL[4],e2xx[4],kBL[4],/*kup[4],e1up[4],x2frac,x2temp,lr, lrfrac,x2,*/ly[4], rest[11], tet, tpt, tst, Bup[4],Bi[4],u[4],uKS[4],uloKS[4],ulo[4],
	kupKS[4],e1upKS[4],BupKS[4],rho,Ttot,sc[4];
doub B[4]={0.,0.,0.,0.}, yyloKS[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
doub knorm,e1norm,ue1,ku,kbpar,kbperp, Be1,Be2,Bn,sin2k,cos2k,Bnu,jVc,aVc,jQc,aQc,jIc,xjIc,xaIc,
aIc,xrVc,rVc,xIc,xQc,xVc,rQc,xrQc,coskb,sinkb,XX,fr,nufr,nW;
int m,k,j,tn,xrn,rn,n,sn,i;

int curr=*(int*)pas;//printf(" %d ", curr);
if((curr>=nthreads)||(curr<0)){
	printf("Error in thread number\n");};
doub nu=ppy[curr].nu;
//t=1.001;
int indx=ppy[curr].indx; doub la=ppy[curr].lamx[0], lb=ppy[curr].lamx[indx]+1e-14;
doub lx,lz=t;int ia=0,ib=indx,ix;
if((lz<=lb) && (la<=lz)){
	while(ib>ia+1){ix=(ia+ib)>>1;lx=ppy[curr].lamx[ix];if(lz<lx){ib=ix;}else{ia=ix;};};}
else{printf("Lookup error 3\n");};
la=ppy[curr].lamx[ia];lb=ppy[curr].lamx[ib];
doub drman=(lz-la)/(lb-la);
if((drman>1.)||(drman<0.)){
//printf("drman error\n");
};
for(m=0;m<4;m++)ly[m]=(1-drman)*ppy[curr].cooxx[m][ia]+drman*ppy[curr].cooxx[m][ib];
for(m=0;m<4;m++)kBL[m]=((1-drman)*ppy[curr].cooxx[m+4][ia]+drman*ppy[curr].cooxx[m+4][ib])/r0;
for(m=0;m<4;m++)e1BL[m]=(1-drman)*ppy[curr].cooxx[m+8][ia]+drman*ppy[curr].cooxx[m+8][ib];

rr=ly[1];costh=ly[2];ph=dphi-ly[3];//metrics is calculated without the exclusion of the polar region!
//remember "-" sign!
if((fabs(rr)>100000.)||(rr<1.)){
printf("Error in geodesic rr=%f\n",rr);};

doub lrz=log(rr);
int indr=ndd-1;
doub lrx,lra=coord[0][0][0],lrb=coord[indr][0][0];ia=0;ib=indr;
if((lrz<=lrb) && (lra<=lrz)){
	while(ib>ia+1){ix=(ia+ib)>>1;lrx=coord[ix][0][0];if(lrz<lrx){ib=ix;}else{ia=ix;};};}
else{printf("Lookup error 4\n");};
lra=coord[ia][0][0];lrb=coord[ib][0][0];
rn=ia;lrman=(lrz-lra)/(lrb-lra);//r coordinates

const doub critan=-((1.-lrman)*theta[rn][0]+lrman*theta[rn+1][0]);//cos(0.07333903);
if(costh>critan){costh=critan;};if(costh<-critan){costh=-critan;};//limiting to the numerical domain
int indth=thlen-1;
doub costhx,tha,thb;ia=0;ib=indth;
while(ib>ia+1){ix=(ia+ib)>>1;costhx=(1.-lrman)*theta[rn][ix]+lrman*theta[rn+1][ix];if(costh<costhx){ib=ix;}else{ia=ix;};};
tha=(1.-lrman)*theta[rn][ia]+lrman*theta[rn+1][ia];thb=(1.-lrman)*theta[rn][ib]+lrman*theta[rn+1][ib];
drman=(costh-tha)/(thb-tha);
tn=thlen-ib-1;thman=1-drman;//theta "traditional" coordinates

/*changed PI to 2*PI*/
phfrac=fmod(ph*phlen/(2*PI)+0.5,phlen);if(phfrac<0)phfrac+=phlen;
k=floor(phfrac);phman=phfrac-k;int ak=(1+k)%phlen;//circular boundary condition 
//if(ph<0){k+=0.;}


if(fdiff==0){sn=0;sing=true;} else {timfrac=(t-mintim)/(maxtim-mintim)*2*fdiff;sn=floor(timfrac);timman=timfrac-sn;sing=false;}
if(t>=maxtim){sn=0;sing=true;};
if(t<=mintim){sn=2*fdiff;sing=true;};


cossq=costh*costh;sinsq=1-cossq;rsq=rr*rr;rhosq=rsq+asq*cossq;Del=rsq-2.*rr+asq;
temp=2.*rr/rhosq;//-,+,+,+ signature covariant KS matrix
iKS[0][0]=temp-1.;iKS[0][1]=temp;iKS[0][2]=0.;iKS[0][3]=-a*temp*sinsq;
iKS[1][0]=iKS[0][1]; iKS[1][1]=1.+temp;iKS[1][2]=0.;iKS[1][3]=-a*(1.+temp)*sinsq;
iKS[2][0]=iKS[0][2];iKS[2][1]=iKS[1][2];iKS[2][2]=rhosq;iKS[2][3]=0.;
iKS[3][0]=iKS[0][3];iKS[3][1]=iKS[1][3];iKS[3][2]=iKS[2][3];iKS[3][3]=sinsq*(rhosq+asq*(1+temp)*sinsq);


if(rr<=rcut){
if(sing)for(m=0;m<11;m++)rest[m]=(1-lrman)*(1-thman)*(1-phman)*(*uu[sn])[k][tn][rn][m]+lrman*(1-thman)*(1-phman)*(*uu[sn])[k][tn][rn+1][m]+(1-lrman)*thman*(1-phman)*(*uu[sn])[k][tn+1][rn][m]+lrman*thman*(1-phman)*
(*uu[sn])[k][tn+1][rn+1][m]+(1-lrman)*(1-thman)*phman*(*uu[sn])[ak][tn][rn][m]+lrman*(1-thman)*phman*(*uu[sn])[ak][tn][rn+1][m]+(1-lrman)*thman*phman*(*uu[sn])[ak][tn+1][rn][m]+lrman*thman*phman*(*uu[sn])[ak][tn+1][rn+1][m];
else for(m=0;m<11;m++)rest[m]=(1-timman)*((1-lrman)*(1-thman)*(1-phman)*(*uu[sn])[k][tn][rn][m]+lrman*(1-thman)*(1-phman)*(*uu[sn])[k][tn][rn+1][m]+(1-lrman)*thman*(1-phman)*(*uu[sn])[k][tn+1][rn][m]+lrman*thman*(1-phman)*
(*uu[sn])[k][tn+1][rn+1][m]+(1-lrman)*(1-thman)*phman*(*uu[sn])[ak][tn][rn][m]+lrman*(1-thman)*phman*(*uu[sn])[ak][tn][rn+1][m]+(1-lrman)*thman*phman*(*uu[sn])[ak][tn+1][rn][m]+lrman*thman*phman*(*uu[sn])[ak][tn+1][rn+1][m])
+timman*((1-lrman)*(1-thman)*(1-phman)*(*uu[sn+1])[k][tn][rn][m]+lrman*(1-thman)*(1-phman)*(*uu[sn+1])[k][tn][rn+1][m]+(1-lrman)*thman*(1-phman)*(*uu[sn+1])[k][tn+1][rn][m]+lrman*thman*(1-phman)*
(*uu[sn+1])[k][tn+1][rn+1][m]+(1-lrman)*(1-thman)*phman*(*uu[sn+1])[ak][tn][rn][m]+lrman*(1-thman)*phman*(*uu[sn+1])[ak][tn][rn+1][m]+(1-lrman)*thman*phman*(*uu[sn+1])[ak][tn+1][rn][m]+lrman*thman*phman*(*uu[sn+1])[ak][tn+1][rn+1][m]);
//{rest[m]=((4.-lrman-thman-phman-timman)*(*uu[sn])[k][tn][rn][m]+lrman*(*uu[sn])[k][tn][rn+1][m]+thman*(*uu[sn])[k][tn+1][rn][m]+phman*(*uu[sn])[ak][tn][rn][m]+timman*(*uu[sn+1])[k][tn][rn][m])/4.;}

//if(rest[1]/rest[0]>1.)rest[1]=rest[0]*1.;//limiting temperature - mostly works in the equatorial plane for high frequency :(
rho=rest[0]*rhonor; Ttot=rest[1]*mp*cc*cc/3/kb/rest[0];
//testing k values
//if((40<k)&&(k<80)){rho/=50.;}//test passed - the sense of rotation is clockwise for 1.57<th<3.14

u[0]=rest[4];u[1]=u[0]*rest[5];u[2]=u[0]*rest[6];u[3]=u[0]*rest[7];
Bi[1]=Bnor*rest[8];Bi[2]=Bnor*rest[9];Bi[3]=Bnor*rest[10];//Btemp=Bnor;

	for(i=0;i<4;i++)for(k=0;k<4;k++){
//	m=26+4*i+k;iMKS[i][k]=(1-lrman)*(1-thman)*(*usgread)[tn][rn][m]+lrman*(1-thman)*(*usgread)[tn][rn+1][m]+(1-lrman)*thman*(*usgread)[tn+1][rn][m]+lrman*thman*(*usgread)[tn+1][rn+1][m];
//	m=47+4*i+k;MKStoKS[i][k]=(1-lrman)*(1-thman)*(*usgread)[tn][rn][m]+lrman*(1-thman)*(*usgread)[tn][rn+1][m]+(1-lrman)*thman*(*usgread)[tn+1][rn][m]+lrman*thman*(*usgread)[tn+1][rn+1][m];
	MKStoKS[i][k]=(1-lrman)*(1-thman)*dxdxp[rn][tn][i][k]+lrman*(1-thman)*dxdxp[rn+1][tn][i][k]+(1-lrman)*thman*dxdxp[rn][tn+1][i][k]+lrman*thman*dxdxp[rn+1][tn+1][i][k];
	}
	for(i=0;i<4;i++)for(k=0;k<4;k++)iMKS[i][k]=0.;
	for(i=0;i<4;i++)for(k=0;k<4;k++)for(j=0;j<4;j++)for(m=0;m<4;m++)iMKS[i][k]+=MKStoKS[m][i]*iKS[m][j]*MKStoKS[j][k];

for(i=1;i<4;i++){uKS[i]=0.;for(k=0;k<4;k++)uKS[i]+=MKStoKS[i][k]*u[k];}

doub expr=fabs((uKS[1]*iKS[0][1] + uKS[3]*iKS[0][3])* (uKS[1]*iKS[0][1] + uKS[3]*iKS[0][3]) - iKS[0][0]*(1. + uKS[1]*uKS[1]*iKS[1][1] + 2.*uKS[1]*uKS[3]*iKS[1][3] + uKS[2]*uKS[2]*iKS[2][2] + uKS[3]*uKS[3]*iKS[3][3]));
uKS[0]=-((uKS[1]*iKS[0][1] + uKS[3]*iKS[0][3] + sqrt(expr))/iKS[0][0]);
for (m=0;m<4;m++){uloKS[m]=0.;for(j=0;j<4;j++)uloKS[m]+=iKS[m][j]*uKS[j];};//uloKS
kupKS[0]=kBL[0]+2*rr*kBL[1]/Del;kupKS[1]=kBL[1];kupKS[2]=kBL[2];kupKS[3]=kBL[3]+a*kBL[1]/Del;//kupKS
e1upKS[0]=e1BL[0]+2*rr*e1BL[1]/Del;e1upKS[1]=e1BL[1];e1upKS[2]=e1BL[2];e1upKS[3]=e1BL[3]+a*e1BL[1]/Del;//e1upKS
ku=0.;ue1=0.;for(m=0;m<4;m++){ue1+=uloKS[m]*e1upKS[m];ku+=uloKS[m]*kupKS[m];};for(m=0;m<4;m++)e1upKS[m]-=ue1/ku*kupKS[m];//update e1upKS

u[0]=uKS[0];//works in simple MKS, where only spatial transformation is performed
for (m=0;m<4;m++){ulo[m]=0.;for(j=0;j<4;j++)ulo[m]+=iMKS[m][j]*u[j];};//ulo

doub yyKS[4][4]={{uKS[0], uKS[1], uKS[2], uKS[3]},{uKS[0]*uloKS[1], -(uKS[0]*uloKS[0] + uKS[3]*uloKS[3]), 0, uKS[3]*uloKS[1]},
{uloKS[2]*uKS[0], uloKS[2]*uKS[1],uloKS[2]*uKS[2]+1, uloKS[2]*uKS[3]},{-uloKS[3]/uloKS[0], 0, 0, 1}},
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
for(m=0;m<4;m++){for(j=0;j<4;j++)B[m]+=yyloKS[m][j]*BupKS[j];};
m=1;
} else {
	uKS[0]=1.;uKS[1]=0.;uKS[2]=0.;uKS[3]=0.;
	doub expr=fabs((uKS[1]*iKS[0][1] + uKS[3]*iKS[0][3])* (uKS[1]*iKS[0][1] + uKS[3]*iKS[0][3]) - iKS[0][0]*(1. + uKS[1]*uKS[1]*iKS[1][1] + 2.*uKS[1]*uKS[3]*iKS[1][3] + uKS[2]*uKS[2]*iKS[2][2] + uKS[3]*uKS[3]*iKS[3][3]));
	uKS[0]=-((uKS[1]*iKS[0][1] + uKS[3]*iKS[0][3] + sqrt(expr))/iKS[0][0]);
	for (m=0;m<4;m++){uloKS[m]=0.;for(j=0;j<4;j++)uloKS[m]+=iKS[m][j]*uKS[j];};//uloKS
	kupKS[0]=kBL[0]+2*rr*kBL[1]/Del;kupKS[1]=kBL[1];kupKS[2]=kBL[2];kupKS[3]=kBL[3]+a*kBL[1]/Del;//kupKS
	e1upKS[0]=e1BL[0]+2*rr*e1BL[1]/Del;e1upKS[1]=e1BL[1];e1upKS[2]=e1BL[2];e1upKS[3]=e1BL[3]+a*e1BL[1]/Del;//e1upKS
	ku=0.;ue1=0.;for(m=0;m<4;m++){ue1+=uloKS[m]*e1upKS[m];ku+=uloKS[m]*kupKS[m];};for(m=0;m<4;m++)e1upKS[m]-=ue1/ku*kupKS[m];//update e1upKS

	u[0]=uKS[0];//works in simple MKS, where only spatial transformation is performed
	for (m=0;m<4;m++){ulo[m]=0.;for(j=0;j<4;j++)ulo[m]+=iMKS[m][j]*u[j];};//ulo

	doub yyKS[4][4]={{uKS[0], uKS[1], uKS[2], uKS[3]},{uKS[0]*uloKS[1], -(uKS[0]*uloKS[0] + uKS[3]*uloKS[3]), 0, uKS[3]*uloKS[1]},
	{uloKS[2]*uKS[0], uloKS[2]*uKS[1],uloKS[2]*uKS[2]+1, uloKS[2]*uKS[3]},{-uloKS[3]/uloKS[0], 0, 0, 1}},
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
for(m=0;m<5;m++)rest[m]=(1-thman)*(1-phman)*uext[k][tn][m]+thman*(1-phman)*uext[k][tn+1][m]+(1-thman)*phman*uext[ak][tn][m]+thman*phman*uext[ak][tn+1][m];
rho=rest[3]*pow(rr/rcut,-rhopo);Ttot=rest[4]*pow(rr/rcut,-1.0);
temp=pow(rr/rcut,-0.5*(1.+rhopo));B[1]=rest[0]*temp;B[2]=rest[1]*temp;B[3]=rest[2]*temp;
}

if(1-fabs(costh)<thlimit){Ttot*=0.000001;};//limiting the numerical domain - works really well, but unphysical!

if(Ttot>maxT){Ttot=maxT;};
if(Ttot<minT){Ttot=minT;};
int indT=stNt; doub Ta=ts[0], Tb=ts[indT];
doub Tx,Tz=Ttot;ia=0;ib=indT;ix;
if((Ta<=Tz) && (Tz<=Tb)){
	while(ib>ia+1){ix=(ia+ib)>>1;Tx=ts[ix];if(Tz<Tx){ib=ix;}else{ia=ix;};};}
else{
	printf("Lookup error 2\n");};
Ta=ts[ia];Tb=ts[ib];
drman=(Tz-Ta)/(Tb-Ta);
//tst=(1-drman)*ts[ia]+drman*ts[ib];
tpt=(1-drman)*tp[ia]+drman*tp[ib];
tet=(1-drman)*te[ia]+drman*te[ib];

if(tet<0){
	tet=0.0001; printf("te<0");};

	doub kxx[4]={0.,0.,0.,0.};for(m=0;m<4;m++){for(j=0;j<4;j++)kxx[m]+=yyloKS[m][j]*kupKS[j];};
	doub e1xx[4]={0.,0.,0.,0.};for(m=0;m<4;m++){for(j=0;j<4;j++)e1xx[m]+=yyloKS[m][j]*e1upKS[j];};
	doub test[4]={0.,0.,0.,0.};for(m=0;m<4;m++){for(j=0;j<4;j++)test[m]+=yyloKS[m][j]*uKS[j];};

knorm=0.;e1norm=0.;for(m=1;m<4;m++){knorm+=kxx[m]*kxx[m];e1norm+=e1xx[m]*e1xx[m];}//kxx[0]=0;//orthonormalization of k with u vector
knorm=-sqrt(knorm);e1norm=sqrt(e1norm);for(m=1;m<4;m++){kxx[m]/=knorm;e1xx[m]/=e1norm;};
e2xx[0]=0.;e2xx[1]=-e1xx[2]*kxx[3]+e1xx[3]*kxx[2];e2xx[2]=-e1xx[3]*kxx[1] + e1xx[1]*kxx[3];e2xx[3]=e1xx[2]*kxx[1] - e1xx[1]*kxx[2];//e2=k x e1 - opposite to Melrose

kbpar=(kxx[1]*B[1]+kxx[2]*B[2]+kxx[3]*B[3]);
kbperp=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]-kbpar*kbpar);

fr=kxx[0];
if(fr<0){
	fr=-fr;//printf("freq=%f ",fr);
};//guard against errors
//tet=8e10/sqrt(rr);//for tests

double magn=(B[1]*B[1]+B[2]*B[2]+B[3]*B[3])/4/PI/mp/rho/cc/cc;
//if(magn>10)tet=1.;//excluding the region of very high magnetization by effectively setting Te=0; magn>30
if(isBcut){if(((rr<9) && (magn+20./(9-rg)*(rr-rg))>30)||((rr>=9)&& (magn>10.))){tet=1.;};};
if(isBred){tet*=exp(-magn/10.);};
//for thickdisk7: magn>30 at r=rg; magn>10 at r=9 => magn-20/(9-rg)*(r-rg)>30

nufr=nu*fr;nW=2*me*cc*2*PI*nufr/(3*ee*kbperp); Bnu=(2*kb*nufr*nufr*tet)/cc/cc;
#include "emis.cpp"
//printf("%.10f ",kbpar);
jIc=(sqrt(3.)*ee*ee*ee*rho*rgrav/2.)/(4.*PI*cc*cc*me)*kbperp*xIc;//unit distance is M=rgrav/2 - this emissivity is already per unit distance
jQc=(sqrt(3.)*ee*ee*ee*rho*rgrav/2.)/(4.*PI*cc*cc*me)*kbperp*xQc;
jVc=fljVc*(ee*ee*ee*rho*rgrav/2.)/(sqrt(3.)*PI*cc*cc*me)*kbpar*xVc;//plus sign, because k vector is along geodesic from the observer to BH, and rad.transfer goes the other way
Bn=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);coskb=kbpar/Bn;sinkb=sqrt(1-coskb*coskb);
XX=tet/The*sqrt(sqrt(2.)*kbperp*1000.*ee/(2.*cc*me*nufr*PI));
rQc=flrQc*kbperp*kbperp*rgrav/2.*rho*ee*ee*ee*ee/cc/cc/cc/me/me/me/4/PI/PI/nufr/nufr/nufr*xrQc*(2.011*exp(-pow(XX,(doub)1.035)/4.7)-cos(XX/2.)*exp(-pow(XX,(doub)1.2)/2.73)-0.011*exp(-XX/47.2));
rVc=flrVc*ee*ee*ee*rho*rgrav/2.*kbpar/PI/cc/cc/me/me/nufr/nufr*xrVc*(1.-0.11*log(1.+0.035*XX));
aIc=jIc/Bnu;aQc=jQc/Bnu;aVc=jVc/Bnu;
Be1=B[1]*e1xx[1]+B[2]*e1xx[2]+B[3]*e1xx[3];Be2=B[1]*e2xx[1]+B[2]*e2xx[2]+B[3]*e2xx[3];
sin2k=-2.*Be1*Be2/(Be1*Be1+Be2*Be2);cos2k=(Be2*Be2-Be1*Be1)/(Be1*Be1+Be2*Be2);

doub The=kb*tet/me/cc/cc,nuc=ee*kbperp/me/cc/2./PI,nus=2./9.*nuc*The*The,//warning: sin(theta) moved from nus to nuc
Xe=nu/nus; xjIc=rho*rgrav/2.*sqrt(2.)*PI*ee*ee*nus/6./The/The/cc*Xe*exp(-pow(Xe,(doub)0.3333));
doub tem=(sqrt(Xe)+pow(2.,11./12.)*pow(Xe,(doub)0.166666));
xjIc=rho*rgrav/2.*sqrt(2.)*PI*ee*ee*nus/6./The/The/cc*tem*tem*exp(-pow(Xe,(doub)0.3333));
xaIc=xjIc/Bnu;
if(xjIc>1e-2){ 
	m=1;}
if(rr<5.){ 
	m=1;}

/*else {
jIc=1.e-20;
xjIc=1.e-30;fr=1.;
jQc=0.;
jVc=0.;
Bn=1.;coskb=0.5;sinkb=sqrt(1-coskb*coskb);
XX=0.;
rQc=0.;
rVc=0.;Bnu=1.e+10;
aIc=jIc/Bnu;aQc=jQc/Bnu;aVc=jVc/Bnu;
xaIc=xjIc/Bnu;
sin2k=0.5;cos2k=sqrt(1-sin2k*sin2k);
};
*/

if(jIc!=jIc){
printf("Epic fail\n");}