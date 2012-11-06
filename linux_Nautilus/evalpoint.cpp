//using namespace std;
bool ju=true;
doub rr, costh, lr, x2,x2frac,x2man, lrfrac,lrman, temp,tempsq, cossq,sinsq,rsq,Del,rhosq,sisq;
doub yylo[4][4],g[4][4], ig[4][4], e1up[4],e1do[4],e2xx[4],kup[4],ly[4], rest[18], tet, tpt, B[4],Bsq[4],Brms[4],u[4],ulo[4],rho,sc[4];
doub kdo[4], knorm,e1norm,ke1,ue1,ku,kbpar,kbperp,e2norm, ke2,e1e2,ll,Bl,Fl,Jl,Hl, Be1,Be2,Bn,sin2k,cos2k,Bnu,jVc,aVc,jQc,aQc,jIc,
aIc,xrVc,rVc,xIc,xQc,xVc,rQc,xrQc,coskb,sinkb,XX,fr,nufr,nW;
int m,k,j,cou;poin *pp; pp=(poin*)pas;doub nu=(*pp).nu;

for(m=0;m<4;m++)ly[m] = gsl_spline_eval ((*pp).xk[m], t, (*pp).acxk[m]);
for(m=0;m<4;m++)kup[m] = gsl_spline_eval ((*pp).xk[m+4], t, (*pp).acxk[m+4])/r0;
for(m=0;m<4;m++)e1up[m] = gsl_spline_eval ((*pp).xk[m+8], t, (*pp).acxk[m+8]);

//if(kup[0]<0){printf("t=%.8f; tmin=%.7f, %Xh, N=%d, kup0(t)=%.3e\n",llog,(*pp).llmax,pp,omp_get_thread_num(),gsl_spline_eval ((*pp).xk[4], t-0.0001, (*pp).acxk[4])/r0);};

//lookup the position, tangential & perpendicular vectors
rr=ly[1];costh=ly[2];
if(costh>cos(0.07333903)){costh=cos(0.07333903);};if(costh<cos(3.06825362)){costh=cos(3.06825362);};//limiting to the numerical domain

lr=log(rr-off); x2 = gsl_spline_eval (th2x2, costh, (*pp).acc);
x2frac=(x2-x2min)/(x2max-x2min)*(thlen-1.);j=floor(x2frac);x2man=x2frac-j;
//if(rr>1e10){printf("rr=%.2e; t=%.2e; N=%d\n",rr,t,omp_get_thread_num());};

if(rr<=rcut){lrfrac=(lr-lrmin)/(lrmax-lrmin)*(rlen-1.);k=floor(lrfrac);lrman=lrfrac-k;
temp=(1. + lrman*(-1. + x2man) - x2man);
for(m=0;m<=17;m++)
{rest[m]=temp*xx[k][j][m] + x2man*xx[k][1 + j][m]
	+ lrman*(xx[1 + k][j][m] + x2man*(-xx[k][1 + j][m] - xx[1 + k][j][m]+ xx[1 + k][1 + j][m]));
};rho=rest[2];rho*=rhonor; rest[3]*=mp*cc*cc/3/kb/rest[2];//correct temperature! was rest[3]*=rhonor/rhocon*mp*cc*cc/3/kb;
u[0]=rest[4];u[1]=rest[5];u[2]=rest[6];u[3]=rest[7];
B[1]=rest[12];B[2]=rest[13];B[3]=rest[14];//linear magnetic fields
Bsq[1]=rest[15];Bsq[2]=rest[16];Bsq[3]=rest[17];//squared magnetic fields
B[1]*=Bnor;B[2]*=Bnor;B[3]*=Bnor;
Bsq[1]*=Bnor*Bnor;Bsq[2]*=Bnor*Bnor;Bsq[3]*=Bnor*Bnor;
}
else {//printf("Second case; r > %f \n",rcut);
	for(m=1;m<=17;m++)rest[m]=(1. - x2man)*xx[ncut-1][j][m] + x2man*xx[ncut-1][1 + j][m];
	rest[0]=rr;rho=rest[2];
	rho*=rhonor*pow(rr/rcut,-rhopo);rest[3]*=mp*cc*cc/3/kb/rest[2]*pow(rr/rcut,-Upo);//correct temperature
	u[0]=rest[4];u[1]=rest[5];u[2]=rest[6];u[3]=rest[7];
	temp=sqrt(rcut/rr);tempsq=rcut/rr;u[0]=1+(u[0]-1)*tempsq;u[1]*=temp;u[2]*=temp*tempsq;u[3]*=temp*tempsq;
//{u[0] -> 1 + a/r, u[1] -> a/Sqrt[r], u[2] -> a/r^(3/2), u[3] -> a/r^(3/2)}
	temp=Bnor*pow(rr/rcut,-Bpo);tempsq=temp*temp;
	B[1]=rest[12];B[2]=rest[13];B[3]=rest[14];//linear magnetic fields
	Bsq[1]=rest[15];Bsq[2]=rest[16];Bsq[3]=rest[17];//squared magnetic fields
	B[1]*=temp;B[2]*=temp;B[3]*=temp;
	Bsq[1]*=tempsq;Bsq[2]*=tempsq;Bsq[3]*=tempsq;};
for(m=1;m<=3;m++)Brms[m]=_copysign(sqrt(Bsq[m]),B[m]);
tet=gsl_spline_eval (tespl, rest[3], (*pp).teacc); tpt = gsl_spline_eval (tpspl, rest[3], (*pp).tpacc);

//u[0]=1.02;u[1]=0.001;u[2]=0.00123;u[3]=0.0006966915908847646;rr=101.860637664794922, costh=cos(1.1711983680725098); asq=0.;a=0.;

cossq=costh*costh;sinsq=1-cossq;rsq=rr*rr;
Del=rsq-2.*rr+asq; rhosq=rsq+asq*cossq; sisq=(rsq+asq)*(rsq+asq)-asq*Del*sinsq;
/*g[0][0]=-sisq/rhosq/Del;//-,+,+,+ signature//contravariant matrix g^munu
g[0][1]=0.;g[0][2]=0.;g[0][3]=-2.*a*rr/rhosq/Del;
g[1][0]=0.;g[1][2]=0.;g[1][3]=0.;g[1][1]=Del/rhosq;g[2][0]=0.;g[2][1]=0.;g[2][3]=0.;
g[2][2]=1./rhosq;g[3][1]=0.;g[3][2]=0.;g[3][0]=g[0][3];g[3][3]=(Del-asq*sinsq)/rhosq/Del/sinsq;
*/

ig[0][0]=2.*rr/rhosq-1.;//-,+,+,+ signature//covariant matrix g_munu
ig[0][1]=0.;ig[0][2]=0.;ig[0][3]=-2.*a*rr*sinsq/rhosq;
ig[1][0]=0.;ig[1][2]=0.;ig[1][3]=0.;ig[1][1]=rhosq/Del;ig[2][0]=0.;ig[2][1]=0.;
ig[2][3]=0.;ig[2][2]=rhosq;ig[3][1]=0.;ig[3][2]=0.;ig[3][0]=ig[0][3];ig[3][3]=sinsq*sisq/rhosq;

u[0]=-((ig[0][3]*u[3] + sqrt(ig[0][3]*ig[0][3]*u[3]*u[3] - ig[0][0]*(1 +
ig[1][1]*u[1]*u[1] + ig[2][2]*u[2]*u[2] + ig[3][3]*u[3]*u[3]))))/ig[0][0];//reconstruction of u[0]

ulo[0]=ig[0][0]*u[0] + ig[0][3]*u[3];ulo[1]=ig[1][1]*u[1];ulo[2]=ig[2][2]*u[2];ulo[3]=ig[0][3]*u[0] + ig[3][3]*u[3];
kdo[0]=ig[0][0]*kup[0] + ig[0][3]*kup[3];kdo[1]=ig[1][1]*kup[1];kdo[2]=ig[2][2]*kup[2];kdo[3]=ig[0][3]*kup[0] + ig[3][3]*kup[3];
ku=0.;ue1=0.;for(m=0;m<4;m++){ue1+=ulo[m]*e1up[m];ku+=ulo[m]*kup[m];};
for(m=0;m<4;m++)e1up[m]-=ue1/ku*kup[m];//update e1up, calculate e1do
e1do[0]=ig[0][0]*e1up[0] + ig[0][3]*e1up[3];e1do[1]=ig[1][1]*e1up[1];e1do[2]=ig[2][2]*e1up[2];e1do[3]=ig[0][3]*e1up[0] + ig[3][3]*e1up[3];

knorm=0.;e1norm=0.;ke1=0.;ku=0.;ue1=0.;
for(m=0;m<4;m++){knorm+=kdo[m]*kup[m]; e1norm+=e1do[m]*e1up[m]; ke1+=e1do[m]*kup[m];ue1+=ulo[m]*e1up[m];ku+=ulo[m]*kup[m];};
//printf("r= %e, norm(k) = %e, norm(e1) = %e, k.e1 = %e, u.e1 = %e, k.u = %e\n",r, knorm,e1norm,ke1,ue1,ku);

/*ll = ulo[3]/ulo[0];Fl=-((ll*ig[0][0] - ig[0][3])/(ll*ig[0][3] - ig[3][3]));
Bl = -((ulo[0] + Fl*ulo[3])/ulo[1]);Hl=(-ig[0][0] - 2*Fl*ig[0][3] - Fl*Fl*ig[3][3])/(Bl*ig[1][1]);
Jl=(-ulo[0] - Hl*ulo[1] - Fl*ulo[3])/ulo[2];doub yy[4][4]={{u[0],u[1],u[2],u[3]},{-1,-Bl,0,-Fl},{1,Hl,Jl,Fl},{-ll,0,0,1}};//changed sign of yy[1]
*/

doub yy[4][4]={{u[0], u[1], u[2], u[3]},{-u[0]*ulo[1]*ig[0][0], (u[0]*ulo[0] + u[3]*ulo[3])/ig[1][1], 0, -u[3]*ulo[1]*ig[0][0]},
{ulo[2]*u[0], ulo[2]*u[1],ulo[2]*u[2]+1, ulo[2]*u[3]},{-ulo[3]/ulo[0], 0, 0, 1}};

for(m=0;m<4;m++){yylo[m][0]=ig[0][0]*yy[m][0] + ig[0][3]*yy[m][3];yylo[m][1]=ig[1][1]*yy[m][1];yylo[m][2]=ig[2][2]*yy[m][2];
yylo[m][3]=ig[0][3]*yy[m][0] + ig[3][3]*yy[m][3];};
for (m=0;m<4;m++){sc[m]=0.;for(j=0;j<4;j++)sc[m]+=yylo[m][j]*yy[m][j];
sc[m]=sqrt(abs(sc[m]));for(j=0;j<4;j++){yylo[m][j]/=sc[m];yy[m][j]/=sc[m];};};
for(j=0;j<4;j++)yylo[0][j]=-yylo[0][j];

/*xx[0] = {u[0], u[1], u[2], u[3]};xx[1] = {1, B, 0, F};xx[2] = {1, H, J, F};xx[3] = {-ll, 0, 0, 1};*/
doub kxx[4]={0.,0.,0.,0.};for(m=0;m<4;m++){for(j=0;j<4;j++)kxx[m]+=yylo[m][j]*kup[j];};
doub e1xx[4]={0.,0.,0.,0.};for(m=0;m<4;m++){for(j=0;j<4;j++)e1xx[m]+=yylo[m][j]*e1up[j];};
doub test[4]={0.,0.,0.,0.};for(m=0;m<4;m++){for(j=0;j<4;j++)test[m]+=yylo[m][j]*u[j];};

knorm=0.;for(m=1;m<4;m++)knorm+=kxx[m]*kxx[m];knorm=sqrt(knorm);//kxx[0]=0;//orthonormalization of k with u vector
for(m=1;m<4;m++){kxx[m]/=knorm;};
e2xx[0]=0.;e2xx[1]=e1xx[2]*kxx[3]-e1xx[3]*kxx[2];e2xx[2]=e1xx[3]*kxx[1] - e1xx[1]*kxx[3];e2xx[3]=-e1xx[2]*kxx[1] + e1xx[1]*kxx[2];//e2=e1 x k - vector product(Melrose p.182)

e2norm=0.;e1norm=0.;knorm=0.;ke1=0.;ke2=0.;e1e2=0.;
//for(m=1;m<4;m++){knorm+=kxx[m]*kxx[m];e1norm+=e1xx[m]*e1xx[m];ke1+=kxx[m]*e1xx[m];ke2+=kxx[m]*e2xx[m];e2norm+=e2xx[m]*e2xx[m];e1e2+=e1xx[m]*e2xx[m];};

//printf("r= %e, norm(k) = %e, norm(e1) = %e, norm(e2) = %e, k.e1 = %e, k.e2 = %e,  e1.e2 = %e\n", r, knorm,e1norm, e2norm,ke1, ke2, e1e2);//k and e1 in the local rest frame: orthonormalized

kbpar=kxx[1]*Brms[1]+kxx[2]*Brms[2]+kxx[3]*Brms[3];
kbperp=sqrt(Bsq[1]+Bsq[2]+Bsq[3]-kbpar*kbpar);
fr=kxx[0];nufr=nu*fr;nW=2*me*cc*2*PI*nufr/(3*ee*kbperp); Bnu=(2*kb*nufr*nufr*tet)/cc/cc;
//if(fr<0){printf("kxx=%.3e; ll=%.7f;t=%.8f; N=%d, kup0(t)=%.3e\n",kxx[0],llog,t,omp_get_thread_num(),gsl_spline_eval ((*pp).xk[4], t, (*pp).acxk[4])/r0);};
#include "emis.cpp"
jIc=(sqrt(3.)*ee*ee*ee*rho*rgrav/2.)/(4*cc*cc*me*PI)*kbperp*xIc;//unit distance is M=rgrav/2 - this emissivity is already per unit distance
jQc=(sqrt(3.)*ee*ee*ee*rho*rgrav/2.)/(4*cc*cc*me*PI)*kbperp*xQc;
jVc=-(ee*ee*ee*rho*rgrav/2.)/(sqrt(3.)*cc*cc*me*PI)*kbpar*xVc;//minus sign, because k vector is along geodesic from the observer to BH, and rad.transfer goes the other way
Bn=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);coskb=(kxx[1]*B[1]+kxx[2]*B[2]+kxx[3]*B[3])/Bn;sinkb=sqrt(1-coskb*coskb);
XX=tet/The*sqrt(sqrt(2.)*sinkb*1000.*Bn*ee/(2.*cc*me*nufr*PI));
rQc=-Bn*Bn*sinkb*sinkb*rgrav/2.*rho*ee*ee*ee*ee/cc/cc/cc/me/me/me/4/PI/PI/nufr/nufr/nufr*xrQc*
(2.011*exp(-pow(XX,(doub)1.035)/4.7)-cos(XX/2.)*exp(-pow(XX,(doub)1.2)/2.73)-0.011*exp(-XX/47.2));
rVc=Bn*ee*ee*ee*rho*rgrav/2.*coskb/PI/cc/cc/me/me/nufr/nufr*xrVc*(1.-0.11*log(1.+0.035*XX));
aIc=jIc/Bnu;aQc=jQc/Bnu;aVc=jVc/Bnu;
Be1=B[1]*e1xx[1]+B[2]*e1xx[2]+B[3]*e1xx[3];Be2=B[1]*e2xx[1]+B[2]*e2xx[2]+B[3]*e2xx[3];
cos2k=(Be1*Be1-Be2*Be2)/(Be1*Be1+Be2*Be2);sin2k=2*Be1*Be2/(Be1*Be1+Be2*Be2);

//15.39 vs 3.79 for 10  => 1.16s  for "cossq -> end" part
//7.98  vs 3.79 for 10  => 0.42s  for "if(rr<)->lookup of tet" part
//18.16 vs 3.79 for 10  => 1.44s  for 12 gsl_spline_eval's!!!
//34.4 vs 3.8   for 10  => 3.06s  for all evalpoint
//evaluation in this point
//printf("lookup at l=%f: t=%f, r=%f, costh=%f, phi=%f\n",t,y[0]/r0,y[1],y[2],y[3]);
//printf("k   : kt=%f, kr=%f, kmu=%f, kphi=%f\n",kup[0],kup[1],kup[2],kup[3]);
//printf("e1up: et=%f, er=%f, emu=%f, ephi=%f\n",e1up[0],e1up[1],e1up[2],e1up[3]);
//printf("r= %.2f, norm(k) = %.2e, norm(e1) = %f, k.e1 = %.2e\n",rr, knorm,e1norm,ke1);
//if(rl[stN]/rg-1.<1e-2)printf("Falls onto the BH\n"); else printf("Does not fall onto the BH\n");*/