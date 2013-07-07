//#include "rkf45.c"
//#include "rk4.c"

#include "stdafx.h"
#include <stdio.h>
#include <iomanip>
#include <string>
#include <math.h>
#include <ctime>
#include <omp.h>
//#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "step.c"//ode
#include "error.c"
#include "cstd.c"
#include "control.c"
#include "evolve.c"
#include "rk2.c"
#include "infnan.c"
#include "fdiv.c"
#include "view.c"
using namespace std;
#include "win_lin_Jon.c"
//#include "win_lin_ADAFs.c"
//#include "win_lin_new.c"
//#include "win_lin.c"

int geon=0;
const int ndd=650, sflen=13/*10*/,flen=3, thn=50, xylen=256*64, dd=3,wdd=11, maxfield=200,maxco=3000,/*1000 for intensity*/maxst=5000, temlen=301, nWlen=120, Tlen=100, nxy=501, snxy=585;
const doub PI = 4.0*atan(1.0), st=0.05,  maxx=12.1, minn=1.2, step=1e-2, sstep=-0.09, Iint=1e-10, Iang=1.;
const doub nWmin=12000.*pow(1.1,-nWlen/2.), lnWmin=log(nWmin), nWmax=12000.*pow(1.1,nWlen/2),lnWmax=log(nWmax), rrmax=3.4e+5, rhoout=130., Tout=1.5e+7,
Tminr=0.4*pow(1.05,-Tlen),Tmin=0.4,Tmax=0.4*pow(1.05,Tlen),lTminr=log(Tminr),lTmin=log(Tmin),lTmax=log(Tmax),mp=1.67e-24, me=9.1e-28, cc=3.e+10, kb=1.38e-16,
rgrav=1.33e+12, ee=4.8e-10, The=me*cc*cc/kb, year=86400.*365.25,Msun=2.00e+33,
//sftabx[sflen+2][2]={{8.45, 120}, {14.90 , 73.}, {22.50,63.}, {43.00, 46.}, {87.73 , 33.}, {102.00, 24.}, {145.00 , 15.}, {230.86, 13.}, {349.00 , 11.}, {674.00, 8.}, {857.00, 7.},{220.,13.}},
sftab[sflen+5][2]={{8.45, 120}, {14.90 , 73.}, {22.50,63.}, {43.00, 46.}, {87.73, 25.9}, {102., 22.3}, {145., 16.4}, {230.86, 12.2}, {349., 10.3}, {674., 8.8}, {857., 8.6}, {1500.,8.6}, {3000.,8.6},{5000.,8.6},
{10e+3,8.6},{20e+3,8.6},{50e+3,8.6}, {220.,13.}}, freqtab[flen+1][2]={{43.,45.},{87.,30.},{230.86,/*22.2*/12.2},{690.,7.}},
tofit[sflen+1][5]={{8.450, 0.683, 0, 0, -0.2500}, {14.90, 0.871, 0, 0, -0.6200}, {22.50, 0.979, 0.1900, 131.0, 0}, {43.00, 1.135, 0.5500, 94.25, 0},
{87.73, 1.841, 1.420, -4., 0}, {102.0, 1.908, 0, 0, 0},{145.0, 2.275, 0, 0, 0}, {230.9, 2.637, 7.398, 111.5, -1.200}, {349.0, 3.181, 6.499, 146.9, -1.500},
{674.0, 3.286,0, 0, 0}, {857.0, 2.867, 0, 0, 0},{1500.,1.,0.,0.,0.},{3000.,1.,0.,0.,0.},{5000.,1.,0.,0.,0.}};

const doub dFnu[sflen+1]={0.031,0.012,0.015,0.026,0.080, 0.1517, 0.2644, 0.1414, 0.1205, 0.3508, 0.2404}, dCP=0.30; doub dLP[3]={0.50/*100.*/,0.658,0.605}; bool isLP87=true;
doub dEVPA[3]={11.,5.4,2.21};//dFnu=0.33, dCP=0.27, dLP=0.26
const gsl_odeiv_step_type * T, * TT, *Tz;//auxilliary for diff. equations

bool inited=false, iswrite,fl=true, echeck1=false, echeck2=false, echeck3=false, Bdir=true, isBcut=false, isBred=false;//LP%, EVPA, CP%
string fdir,outstr,instr;
int fnum, cas,th_id,i,/*j,*/k, w, sp, ncut, nm,co,mco,fdiff=0,mintim,maxtim,stNt,loaded[maxfield];string fif="";

std::clock_t start;
doub rg, temp, rat, dof=0.,t, fact=1.,xisq, ans, th, tth,xth, r, r0, a, off, res, asq, accur,accurr,maxy,rate,ratex,maxT,minT, xxisq,xang, ss=1e-2, dense=1.,
		fljVc=1., flrQc=1.,flrVc=1.;
doub resid[15][4],delta[15],Jac[15][3],bb[3],matr[3][3],det;
doub usp[maxfield][phlen][thlen][4], uspKS[maxfield][phlen][thlen];
typedef doub (*ausar)[snxy+1][snxy+1][sflen+1][5];ausar ausin = (ausar) new doub[snxy+1][snxy+1][sflen+1][5];ausar ausin2 = (ausar) new doub[snxy+1][snxy+1][5];
typedef double (*arra)[nxy+1][nxy+1][5]; arra intab = (arra) new double[nxy+1][nxy+1][5];
typedef double (*para)[20];para params = (para) new double[20];
doub rtab[2000],Tstab[2000],dphi=0.; //for calculating the electron temperature

float ww[phlen][thlen][rlen][wdd], uext[phlen][thlen][5], dxdxp[ndd][thlen][4][4],coord[ndd][thlen][2],coordx[ndd][thlen][2];//float uu[2*fdiffmax+1][phlen][thlen][rlen][wdd];
//typedef float (*uuarr)[2*fdiffmax+1][phlen][thlen][rlen][wdd];uuarr uu = (uuarr) new float[2*fdiffmax+1][phlen][thlen][rlen][wdd];
typedef float (*uuarr)[phlen][thlen][rlen][wdd];uuarr uu[200];//130 fits in 64GB memory for Jon's simulations from 2012
//typedef float (*uuarr)[phlen][thlen][rlen][wdd]; double *uu[2*fdiffmax+1];
//[2*fdiffmax+1]
typedef double (*usgarr)[thlen][rlen][usgsize];usgarr usgread = (usgarr) new double[thlen][rlen][usgsize];
doub xx[rlen][dd];
doub damp,TpTe,heat,ts[10000], te[10000],tp[10000], rcut, rhopo, Upo, Bpo, rhonor,Bnor,Bcon,Bmax,Bmin,rhocon,rhomax,rhomin,rmin,rmax,Ucon,Umax,Umin,
lrmin,lrmax, thlimit,
theta[ndd][thlen], rad[rlen], lrar[rlen], x2min,x2max, in[sflen+1][5], totin[sflen+1],LPo[sflen+1],CP[sflen+1],
xtotin[sflen+1],xLPo[sflen+1],xCP[sflen+1],xEVPA[sflen+1], ytotin[sflen+1],yLPo[sflen+1],yCP[sflen+1],
ztotin[sflen+1][4],zLPo[sflen+1][4],zCP[sflen+1][4],zEVPA[sflen+1][4],
err[sflen+1],ang[sflen+1],avrate,Te6,avTpTe,avTe6;doub jI[Tlen+1][nWlen+1],jQ[Tlen+1][nWlen+1], jV[Tlen+1][nWlen+1],rQ[2*Tlen+1],rV[2*Tlen+1],cheat,crhonor,lastxisq,
lastrhonor,lastheat;
typedef struct {doub lamx[maxco],cooxx[12][maxco];doub llmin,llmax,nu;int indx;} poinx;
poinx ppy[nthreads];

int geodes (doub t, const doub y[], doub f[], void *params)
{//geodesic equations
#include "geodes.cpp"
}
int solte (doub rz, const doub zz[],doub ff[],void *pas)
{//ff[0]->Te; ff[1]->Tp
int indT=ndd-1; doub irmin=rtab[0], irmax=1.000001*rtab[indT], iTsmin=Tstab[0], iTsmax=Tstab[indT];
doub ra=irmin,rb=irmax,rx;int ia=0,ib=indT,ix;
if((ra<=rz) && (rz<=rb)){
	while(ib>ia+1){ix=(ia+ib)>>1;rx=rtab[ix];if(rz<rx){ib=ix;}else{ia=ix;};};}
else{printf("Lookup error 1\n");};
ra=rtab[ia];rb=rtab[ib];
doub drman=(rz-ra)/(rb-ra);
ts[stNt]=(1-drman)*Tstab[ia]+drman*Tstab[ib];

doub Dts=(Tstab[ib]-Tstab[ia])/(rtab[ib]-rtab[ia]),
dcap=2.-0.49/(0.7+kb*zz[0]/me/cc/cc)/(0.7+kb*zz[0]/me/cc/cc);
doub rhofun=rhocon*pow(rz/rcut,-rhopo),
nuc=8.9e-11*pow(zz[0]/(doub)3e10,(doub)-1.5)*rhofun/1e7,
tail=nuc*2.*mp*rgrav*rgrav*rgrav*rz*rz*rhofun/rate*(zz[1]-zz[0])+
2*(sqrt(zz[1])/(sqrt(zz[1])+heat*sqrt(zz[0]))-heat/dcap*sqrt(zz[0])/(sqrt(zz[1])+heat*sqrt(zz[0])))*Dts;
ff[0]=(2*Dts-tail)/(dcap+1);
ff[1]=tail+(2*Dts-tail)/(dcap+1);
return GSL_SUCCESS;}
int init(int sp, int fmin, int fmax, int sep)
{
#include "init.cpp"
}
int trans (doub llog, const doub yyy[], doub ff[], void *pas)
{//polarized radiative transfer
//#include "trans.cpp"
#include "transnew.cpp"
}
int main(int argc, char* argv[])
{int n1,n2,n3;
for(n1=0;n1<2*fdiff+1;n1++)uu[n1]=(uuarr) new float[phlen][thlen][rlen][wdd];
//char * descr = getenv("LSB_JOBINDEX");
char * descr = getenv("PBS_ARRAYID");
int nP,nPeff,kk,j,ittot=0, niter=0, oo,ix,il,iy, kmin,kmax,kstep,sep; doub t;

doub inp[4][3],dheat,drho,dtheta,ddh,ddr,dth;//sp=0,a=0; 1->0.5; 2->0.7; 3->0.9; 4->0.98
T  = gsl_odeiv_step_rk2;accur =3e-4;
TT = gsl_odeiv_step_rk2;accurr=5e-4;r0=20000.;//for real calculations//wfit[0]=false;wfit[1]=false;wfit[2]=true;//LP%, EVPA, CP%

start = clock();
//a - spin, r0 - where we shoot from in radius, th - cosine of polar angle, b - impact parameter, 
//beta - angle in a picture plane from the direction to the north pole counterclockwise //BH rotates clockwise as viewed from above the north pole
lastheat=heat;lastrhonor=rhonor;sep=50;fact=1.;oo=1;iswrite=true;int sear;kmin=4;kmax=10;

#pragma omp parallel num_threads(nthreads) shared(ittot) private(w)
{w = omp_get_thread_num();printf("Hello from thread %d\n", w);}
if(argc==5){sp=atoi(argv[1])-1;co=atoi(argv[2]);Bpo=atof(argv[3]);mco=floor(Bpo+0.001);sear=atoi(argv[4]);} else return 0;

switch (sear){
case 0: //surfs the parameter space near best fit to flux
	//#include "s_space.cpp"
	#include "m_space.cpp"
		break;
case 1: //quick computation of spectrum for certain sp, heat, rhonor, th
	//#include "s_quick.cpp"
	#include "m_quick.cpp"
		break;
case 2: //convergence studies
	//#include "s_conv.cpp"
	#include "m_conv.cpp"
		break;
case 3: //surf region close to best model
	//#include "s_surf.cpp"
	#include "m_surf.cpp"
		break;
case 4: //image of the emitting region
	//#include "s_imag.cpp"
	#include "m_imag.cpp"
		break;
case 5: //searches for minimum with a steepest descent method
	#include "m_sear.cpp"
		break;
case 6: //averages temperature and density profiles to find Te(Ts) and Tp(Ts) functions
	#include "m_ts.cpp"
		break;}

printf ("Iterations per second = % .18f\n", doub(ittot)/ans);
printf ("Iterations = %d\n ", ittot);
return 0;}
