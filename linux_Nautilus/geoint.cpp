int p,stN;doub tpast,y[12],der[12],txx,lam[maxco],coox[12][maxco];//maxco=5000
gsl_odeiv_step *s; gsl_odeiv_control *c; gsl_odeiv_evolve *e;
//poin pp;
//poinx ppx;
int currth; currth = omp_get_thread_num();
s = gsl_odeiv_step_alloc (T, 12);//abs accuracy 1e-6 and method gsl_odeiv_step_rkf45 are the best to generate 200-400 points on a geodesic
c = gsl_odeiv_control_standard_new(0.0, accur, 1.0, 0.0);//c = gsl_odeiv_control_y_new (0.0, accur);
e = gsl_odeiv_evolve_alloc (12);
gsl_odeiv_system sys = {geodes, NULL, 12, NULL};//prepare to solve ODE of 12 variables

rg=1+sqrt(1.-asq);
doub bsq=b*b, sinbeta=sin(beta), cosbeta=cos(beta),sinbetasq=sinbeta*sinbeta,cosbetasq=1-sinbetasq,
 costh=cos(th),sinthsq=1-costh*costh,sinth=sqrt(sinthsq), costhsq=1-sinthsq, r0sq=r0*r0,rmin2=r0-2.;
/*Q = 0.25*r0sq*((-4.*asq*costhsq*((-4. + r0)*r0 + 4.*costhsq))/(rmin2*rmin2) + 4.*bsq*(cosbetasq + costhsq*sinbetasq));Lz = r0*sinth*(b*sinbeta + (2.*asq*sinth)/rmin2);
En = 0.5*((-1.*(bsq*rmin2 + 2.*r0sq*r0))/r0sq + (asq*sinthsq)/rmin2);Lzsq=Lz*Lz; Ensq=En*En;*/ //invariants, calculated for b<<r0
stN=0;
doub tin = 0., rin = 0.5*bsq/r0 + r0, phiin = (-1.*b*(-1.*r0 +b*cosbeta*costh/sinth)/sinth*sinbeta)/r0sq,
muin = (-1.*bsq*costh + 2.*r0sq*costh - 2.*b*r0*cosbeta*sinth)/2./r0sq,
tprin=r0, rprin = bsq/(2.*r0) - r0, muprin = (-1.*bsq*costh - 1.*b*r0*cosbeta*sinth)/r0sq,
phprin=(b*(r0 - 2.*b*cosbeta*costh/sinth)/sinth*sinbeta)/r0sq;//initial position in the picture plane, calculated for b<<r0
doub tain=0.,rain=-(b*cosbeta)/r0,  muain = (2.*b*r0*cosbeta*costh - (bsq-2.*r0sq+2.*bsq*cosbetasq)*sinth)/(2.*r0sq*r0),
pain=(b*costh/sinth*(r0 - 2.*b*cosbeta*costh/sinth)/sinth*sinbeta)/r0sq*r0;

t=0.; doub t1 = 20., tmax;//initial step h, then it changes adaptively, maximum t1=20 is never reached
y[0] = tin; y[1]= phiin; y[2]=rprin; y[3]=muprin; y[4]=rin; y[5]=muin; y[6]=tprin, y[7]=phprin;
y[8]=tain;y[9]=rain;y[10]=muain;y[11]=pain;//slightly higher CP at low f
//y[8]=0.;y[9]=(b*cosbeta)/r0;y[10]=(bsq*(-1. + cosbetasq*(-1. + 2.*costhsq)) - 2.*(-1. + costhsq)*r0sq)/(2.*r0sq*r0*sinthsq);
//y[11]=-(b*costh*sinbeta*(-2.*b*cosbeta*costh +r0*sinth))/(r0sq*r0*sinthsq*sinth);//initial tangential vector of geodesic

geodes (t, y, der, NULL);tpast=t;lam[0]=0.;
coox[0][stN]=y[0];coox[3][stN]=y[1];coox[1][stN]=y[4];coox[2][stN]=y[5];coox[4][stN]=y[6];coox[7][stN]=y[7];coox[5][stN]=y[2];coox[6][stN]=-y[3]/sqrt(1-y[5]*y[5]);
for(p=8;p<12;p++)coox[p][stN]=y[p];
for(p=0;p<12;p++){ppy[currth].cooxx[p][stN]=coox[p][stN];}
ppy[currth].lamx[stN]=lam[stN];

//copy initial points into array - 4 lines not really needed
//y : t, phi, rp ,mup, r, mu ,tpr,phpr
//f : t',phi',r'',mu'',r',mu',t'',ph''
t1=2.1;doub h=step;
while (t < t1)//main integration cycle
{	stN++;
	int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t1, &h, y);lam[stN]=t;//evolve ODE
	//tl[stN]=y[0];phil[stN]=y[1];rl[stN]=y[4];mul[stN]=y[5];ktl[stN]=y[6];kphil[stN]=y[7];krl[stN]=y[2];kmul[stN]=-y[3]/sqrt(1-y[5]*y[5]);
	//e1[0][stN]=y[8];e1[1][stN]=y[9];e1[2][stN]=y[10];e1[3][stN]=y[11];
coox[0][stN]=y[0];coox[3][stN]=y[1];coox[1][stN]=y[4];
//if((y[4]<50.0)&& (stN>0)){	printf("Error\n");}
coox[2][stN]=y[5];coox[4][stN]=y[6];coox[7][stN]=y[7];coox[5][stN]=y[2];coox[6][stN]=-y[3]/sqrt(1-y[5]*y[5]);
for(p=8;p<12;p++)coox[p][stN]=y[p];
for(p=0;p<12;p++){ppy[currth].cooxx[p][stN]=coox[p][stN];}
ppy[currth].lamx[stN]=lam[stN];
	if(tpast==t){printf("Not monotonic!\n");};//needed to check if "float" calculation works
tpast=t; if(h>step)h=step;
if ((status != GSL_SUCCESS) || (y[4]<(1.+ss)*rg)|| (y[4]>r0) || (stN>maxco)) break;}
//check if falls onto BH, flies out of the domain, or exceeds array length
if(stN>maxco)printf("Too many points requested at t=%e\n",t);
//printf("%d ",stN);
tmax=lam[stN];
txx=r0/3.;ppy[currth].llmin=-log(txx+sqrt(1+txx*txx));
txx=r0*(tmax-1.)/3.;ppy[currth].llmax=log(txx+sqrt(1+txx*txx));
ppy[currth].indx=stN;
//for(p=0;p<12;p++)for(k=0;k<=stN;k++){ppy[currth].cooxx[p][k]=coox[p][k];}

//for(k=0;k<=stN;k++){ppy[currth].lamx[k]=lam[k];};
if(stN<10000){
//printf("Hehe");
};
//tlin rlin philin mulin kt kr kmu kphi es[0..3]
gsl_odeiv_evolve_free (e);gsl_odeiv_control_free (c);
gsl_odeiv_step_free (s);