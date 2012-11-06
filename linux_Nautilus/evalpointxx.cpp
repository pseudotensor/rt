bool ju=true;
doub rr, costh, lr, Bsq2B=1., lrfrac,lrman, x2,x2frac,x2man,ph,phfrac,phman, x2temp,temp,tempsq, cossq,sinsq,rsq,Del,rhosq;
doub ig[4][4], e1BL[4],e1up[4],e2xx[4],kBL[4],kup[4],ly[4], rest[11], tet, tpt, Bup[4],Bi[4],u[4],ulo[4],rho,sc[4];
doub knorm,e1norm,ue1,ku,kbpar,kbperp, Be1,Be2,Bn,sin2k,cos2k,Bnu,jVc,aVc,jQc,aQc,jIc,
aIc,xrVc,rVc,xIc,xQc,xVc,rQc,xrQc,coskb,sinkb,XX,fr,nufr,nW;
int m,k,j,n;poin *pp; pp=(poin*)pas;doub nu=(*pp).nu;

for(m=0;m<4;m++)ly[m] = gsl_spline_eval ((*pp).xk[m], t, (*pp).acxk[m]);
for(m=0;m<4;m++)kBL[m] = gsl_spline_eval ((*pp).xk[m+4], t, (*pp).acxk[m+4])/r0;
for(m=0;m<4;m++)e1BL[m] = gsl_spline_eval ((*pp).xk[m+8], t, (*pp).acxk[m+8]);

rr=ly[1];costh=ly[2];ph=ly[3];//metrics is calculated without the exclusion of the polar region!

//rr=1.80051;costh=cos(0.503047);ph=3.09251;//test

lr=log(rr-off); x2 = gsl_spline_eval (th2x2, costh, (*pp).acc);
cossq=costh*costh;sinsq=1-cossq;rsq=rr*rr;rhosq=rsq+asq*cossq;Del=rsq-2.*rr+asq;
/*g[0][0]=-sisq/rhosq/Del;//-,+,+,+ signature//contravariant matrix g^munu
g[0][1]=0.;g[0][2]=0.;g[0][3]=-2.*a*rr/rhosq/Del;
g[1][0]=0.;g[1][2]=0.;g[1][3]=0.;g[1][1]=Del/rhosq;g[2][0]=0.;g[2][1]=0.;g[2][3]=0.;
g[2][2]=1./rhosq;g[3][1]=0.;g[3][2]=0.;g[3][0]=g[0][3];g[3][3]=(Del-asq*sinsq)/rhosq/Del/sinsq;
*/
ig[0][0]=2.*rr/rhosq-1.;//-,+,+,+ signature//covariant MKS matrix g_munu
ig[0][1]=(2.*rr*(-off + rr))/rhosq;ig[0][2]=0.;ig[0][3]=(-2.*a*PI*rr*sinsq)/rhosq;
ig[1][0]=ig[0][1];ig[1][1]=(-off + rr)*(-off + rr)*(1. + (2.*rr)/rhosq);
ig[1][2]=0.;ig[1][3]=-(a*PI*(-off + rr)*(1 + (2*rr)/rhosq)*sinsq);
x2temp=PI*(h1 + 7.*(1 - h1)*pow((doub)(-1. + 2.*x2),(doub)6.));
ig[2][0]=0.;ig[2][1]=0.;ig[2][2]=rhosq*x2temp*x2temp;ig[2][3]=0.;
ig[3][0]=ig[0][3];ig[3][1]=ig[1][3];ig[3][2]=0.;
ig[3][3]=PI*PI*sinsq*(rhosq + asq*(1. + (2.*rr)/rhosq)*sinsq);

doub critan=0.07333903;
if(costh>cos(critan)){costh=cos(critan);};if(costh<cos(PI-critan)){costh=cos(PI-critan);};//limiting to the numerical domain
x2 = gsl_spline_eval (th2x2, costh, (*pp).acc);

x2frac=(x2-x2min)/(x2max-x2min)*(thlen-1.);j=floor(x2frac);x2man=x2frac-j;//j=floor(x2frac+0.5);
phfrac=fmod(ph*phlen/PI+0.5,phlen);if(phfrac<0)phfrac+=phlen;
k=floor(phfrac);phman=phfrac-k;int ak=(1+k)%phlen;//circular boundary condition 

if(rr<=rcut){lrfrac=(lr-lrmin)/(lrmax-lrmin)*(rlen-1.);n=floor(lrfrac);lrman=lrfrac-n; 
for(m=0;m<11;m++)
{//rest[m]=((2.-lrman-x2man)*ww[k][j][n][m]+lrman*ww[k][j][n+1][m]+x2man*ww[k][j+1][n][m])/3.;//23.1s
rest[m]=((3.-lrman-x2man-phman)*ww[k][j][n][m]+lrman*ww[k][j][n+1][m]+x2man*ww[k][j+1][n][m]+phman*ww[ak][j][n][m])/3.;};

rho=rest[0]*rhonor; rest[1]*=mp*cc*cc/3/kb/rest[0];
u[0]=rest[4];u[1]=u[0]*rest[5];u[2]=u[0]*rest[6];u[3]=u[0]*rest[7];
Bi[1]=Bnor*rest[8];Bi[2]=Bnor*rest[9];Bi[3]=Bnor*rest[10];//Btemp=Bnor;
}
else {n=ncut-1;
for(m=0;m<11;m++){
	//rest[m]=(1.-x2man)*ww[k][j][n][m]+x2man*ww[k][j+1][n][m];
	rest[m]=((2.-x2man-phman)*ww[k][j][n][m]+x2man*ww[k][j+1][n][m]+phman*ww[ak][j][n][m])/2.;};

	rho=rest[0]*rhonor*pow(rr/rcut,-rhopo);
	rest[1]*=mp*cc*cc/3/kb/rest[0]*pow(rr/rcut,-Upo);//correct temperature
u[0]=rest[4];u[1]=u[0]*rest[5];u[2]=u[0]*rest[6];u[3]=u[0]*rest[7];
Bi[1]=rest[8];Bi[2]=rest[9];Bi[3]=rest[10];
temp=sqrt(rcut/rr);tempsq=rcut/rr;//u[0]=1+(u[0]-1)*tempsq;
u[1]*=tempsq;u[2]*=temp*tempsq;u[3]*=temp*tempsq;
Bi[1]*=Bnor*tempsq*tempsq*temp;Bi[3]*=Bnor*tempsq*tempsq*temp;Bi[2]*=Bnor*tempsq*temp;};

if(rest[1]>maxT){rest[1]=maxT;};
tet=gsl_spline_eval (tespl, rest[1], (*pp).teacc); tpt = gsl_spline_eval (tpspl, rest[1], (*pp).tpacc);
if(tet<0){
	tet=0;printf("te<0");};
if(averaged)Bsq2B=rest[2]; else Bsq2B=1.;//boost of <B^2>/<B>^2
//u[1]=-0.296765;u[2]= 0.00242941;u[3]=0.0677352;//test
doub expr=(u[1]*ig[0][1] + u[3]*ig[0][3])*
(u[1]*ig[0][1] + u[3]*ig[0][3]) - ig[0][0]*(1 + u[1]*u[1]*ig[1][1] + 2.*u[1]*u[3]*ig[1][3] + 
u[2]*u[2]*ig[2][2] + u[3]*u[3]*ig[3][3]);
u[0]=-((u[1]*ig[0][1] + u[3]*ig[0][3] + sqrt(expr))/ig[0][0]);


for (m=0;m<4;m++){ulo[m]=0.;for(j=0;j<4;j++)ulo[m]+=ig[m][j]*u[j];};//ulo
kup[0]=kBL[0]+2*rr*kBL[1]/Del;kup[1]=kBL[1]/(rr-off);kup[2]=kBL[2]/x2temp;kup[3]=kBL[3]/PI+a*kBL[1]/Del/PI;//kup
e1up[0]=e1BL[0]+2*rr*e1BL[1]/Del;e1up[1]=e1BL[1]/(rr-off);e1up[2]=e1BL[2]/x2temp;e1up[3]=e1BL[3]/PI+a*e1BL[1]/Del/PI;//e1up
ku=0.;ue1=0.;for(m=0;m<4;m++){ue1+=ulo[m]*e1up[m];ku+=ulo[m]*kup[m];};for(m=0;m<4;m++)e1up[m]-=ue1/ku*kup[m];//update e1up
//knorm=0.;e1norm=0.;ke1=0.;ku=0.;ue1=0.;
//for(m=0;m<4;m++){knorm+=kdo[m]*kup[m]; e1norm+=e1do[m]*e1up[m]; ke1+=e1do[m]*kup[m];ue1+=ulo[m]*e1up[m];ku+=ulo[m]*kup[m];};
//printf("r= %e, norm(k) = %e, norm(e1) = %e, k.e1 = %e, u.e1 = %e, k.u = %e\n",r, knorm,e1norm,ke1,ue1,ku);
//kup[0]=4.5118911870585166;kup[1]=-0.69886108066474151;kup[2]=0.25702487876719615;kup[3]=4.1686907691517288;//test

doub yy[4][4]={{u[0], u[1], u[2], u[3]},{u[0]*ulo[1], -(u[0]*ulo[0] + u[3]*ulo[3]), 0, u[3]*ulo[1]},
{ulo[2]*u[0], ulo[2]*u[1],ulo[2]*u[2]+1, ulo[2]*u[3]},{-ulo[3]/ulo[0], 0, 0, 1}}, yylo[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},
xpr[4]={u[0]*ulo[1], -(u[0]*ulo[0] + u[3]*ulo[3]), 0, u[3]*ulo[1]},xlopr[4]={0,0,0,0};
for(j=0;j<4;j++)for(n=0;n<4;n++)xlopr[j]+=ig[j][n]*xpr[n];
for(m=0;m<4;m++)for(j=0;j<4;j++)for(n=0;n<4;n++)yylo[m][n]+=ig[j][n]*yy[m][j];//already transposed matrix

for (m=0;m<4;m++){sc[m]=0.;for(j=0;j<4;j++)sc[m]+=yylo[m][j]*yy[m][j];
sc[m]=sqrt(fabs(sc[m]));for(j=0;j<4;j++){yylo[m][j]/=sc[m];yy[m][j]/=sc[m];};};
for(j=0;j<4;j++)yylo[0][j]=-yylo[0][j];
for(m=0;m<4;m++){yy[1][m]=xpr[m];for(j=0;j<4;j++)yy[1][m]-=yy[3][m]*xlopr[j]*yy[3][j];};
//updating x1 vector - need to do that in MKS ;

for(j=0;j<4;j++){yylo[1][j]=0.;for(n=0;n<4;n++)yylo[1][j]+=ig[n][j]*yy[1][n];};
sc[1]=0.;for(j=0;j<4;j++)sc[1]+=yylo[1][j]*yy[1][j];sc[1]=sqrt(fabs(sc[1]));
for(j=0;j<4;j++){yylo[1][j]/=sc[1];yy[1][j]/=sc[1];}


doub blo=(ulo[1]*Bi[1]+ulo[2]*Bi[2]+ulo[3]*Bi[3]);
Bup[0]=blo;for(m=1;m<4;m++)Bup[m]=(Bi[m]+blo*u[m])/u[0];
doub kxx[4]={0.,0.,0.,0.};for(m=0;m<4;m++){for(j=0;j<4;j++)kxx[m]+=yylo[m][j]*kup[j];};
doub e1xx[4]={0.,0.,0.,0.};for(m=0;m<4;m++){for(j=0;j<4;j++)e1xx[m]+=yylo[m][j]*e1up[j];};
doub test[4]={0.,0.,0.,0.};for(m=0;m<4;m++){for(j=0;j<4;j++)test[m]+=yylo[m][j]*u[j];};
doub B[4]={0.,0.,0.,0.};for(m=0;m<4;m++){for(j=0;j<4;j++)B[m]+=yylo[m][j]*Bup[j];};
//for(m=0;m<4;m++)B[m]*=Bnor;//later

knorm=0.;e1norm=0.;for(m=1;m<4;m++){knorm+=kxx[m]*kxx[m];e1norm+=e1xx[m]*e1xx[m];}//kxx[0]=0;//orthonormalization of k with u vector
knorm=sqrt(knorm);e1norm=sqrt(e1norm);for(m=1;m<4;m++){kxx[m]/=knorm;e1xx[m]/=e1norm;};
e2xx[0]=0.;e2xx[1]=e1xx[2]*kxx[3]-e1xx[3]*kxx[2];e2xx[2]=e1xx[3]*kxx[1] - e1xx[1]*kxx[3];e2xx[3]=-e1xx[2]*kxx[1] + e1xx[1]*kxx[2];//e2=e1 x k - vector product(Melrose p.182)

//e2norm=0.;e1norm=0.;knorm=0.;ke1=0.;ke2=0.;e1e2=0.;
//for(m=1;m<4;m++){knorm+=kxx[m]*kxx[m];e1norm+=e1xx[m]*e1xx[m];ke1+=kxx[m]*e1xx[m];ke2+=kxx[m]*e2xx[m];e2norm+=e2xx[m]*e2xx[m];e1e2+=e1xx[m]*e2xx[m];};

//printf("r= %e, norm(k) = %e, norm(e1) = %e, norm(e2) = %e, k.e1 = %e, k.e2 = %e,  e1.e2 = %e\n", r, knorm,e1norm, e2norm,ke1, ke2, e1e2);//k and e1 in the local rest frame: orthonormalized

//kbpar=kxx[1]*Brms[1]+kxx[2]*Brms[2]+kxx[3]*Brms[3];
//kbperp=sqrt(Bsq[1]+Bsq[2]+Bsq[3]-kbpar*kbpar);
kbpar=(kxx[1]*B[1]+kxx[2]*B[2]+kxx[3]*B[3]);
kbperp=Bsq2B*sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]-kbpar*kbpar);

fr=kxx[0];
if(fr<0){
	fr=-fr;printf("freq=%f ",fr);};//guard against errors

nufr=nu*fr;nW=2*me*cc*2*PI*nufr/(3*ee*kbperp); Bnu=(2*kb*nufr*nufr*tet)/cc/cc;
#include "emis.cpp"
jIc=(sqrt(3.)*ee*ee*ee*rho*rgrav/2.)/(4*cc*cc*me*PI)*kbperp*xIc;//unit distance is M=rgrav/2 - this emissivity is already per unit distance
jQc=(sqrt(3.)*ee*ee*ee*rho*rgrav/2.)/(4*cc*cc*me*PI)*kbperp*xQc;
jVc=-(ee*ee*ee*rho*rgrav/2.)/(sqrt(3.)*cc*cc*me*PI)*kbpar*xVc;//minus sign, because k vector is along geodesic from the observer to BH, and rad.transfer goes the other way
Bn=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);coskb=kbpar/Bn;sinkb=sqrt(1-coskb*coskb);
XX=tet/The*sqrt(sqrt(2.)*sinkb*1000.*Bn*ee/(2.*cc*me*nufr*PI));
rQc=-Bsq2B*Bsq2B*Bn*Bn*sinkb*sinkb*rgrav/2.*rho*ee*ee*ee*ee/cc/cc/cc/me/me/me/4/PI/PI/nufr/nufr/nufr*xrQc*
(2.011*exp(-pow(XX,(doub)1.035)/4.7)-cos(XX/2.)*exp(-pow(XX,(doub)1.2)/2.73)-0.011*exp(-XX/47.2));
rVc=Bn*ee*ee*ee*rho*rgrav/2.*coskb/PI/cc/cc/me/me/nufr/nufr*xrVc*(1.-0.11*log(1.+0.035*XX));
aIc=jIc/Bnu;aQc=jQc/Bnu;aVc=jVc/Bnu;
Be1=B[1]*e1xx[1]+B[2]*e1xx[2]+B[3]*e1xx[3];Be2=B[1]*e2xx[1]+B[2]*e2xx[2]+B[3]*e2xx[3];
cos2k=(Be1*Be1-Be2*Be2)/(Be1*Be1+Be2*Be2);sin2k=2*Be1*Be2/(Be1*Be1+Be2*Be2);
if(jIc!=jIc){
	m=1;}

// for averaged: Bsq2B factor - in total+linear emissivities and Faraday conversion
// pactially factor Bsq2B in circular emissivity, no factor Bsq2B in Faraday rotation

//15.39 vs 3.79 for 10  => 1.16s  for "cossq -> end" part
//7.98  vs 3.79 for 10  => 0.42s  for "if(rr<)->lookup of tet" part
//18.16 vs 3.79 for 10  => 1.44s  for 12 gsl_spline_eval's!!!
//34.4 vs 3.8   for 10  => 3.06s  for all evalpoint