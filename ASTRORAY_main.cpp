/***********************************************************************************
#    Copyright 2014 Roman Shcherbakov
#                   
#              ASTRORAY version 1.0   (released June 25, 2014)
#
#    ASTRORAY v1.0 is a program that performs general relativistic polarized radiative transfer 
#    of synchrotron radiation near black holes. The code employs ray tracing technique 
#    and can easily handle large optical depth (~10,000). 
#    ASTRORAY produces images as well as spectra.
#    
#    The latest version of ASTRORAY code can be retrieved from 
#    our distribution website: http://astroman.org/code/ASTRORAY/.
#
#    This version of ASTRORAY is configured to use input files from 3D version of HARM code,
#    which is maintained by Jonathan McKinney (UMD).
#     The code assumes that the source is a plasma near a
#    black hole described by Kerr-Schild coordinates.
#    Plasma radiates via thermal synchrotron.
#        
#    Please, cite the following paper in any scientific literature 
#    that results from use of any part of ASTRORAY:
# 
#    Shcherbakov R. V., Penna R. F., & McKinney J. C. 2012, Astrophysical Journal, 755, 133
#              http://adsabs.harvard.edu/abs/2012ApJ...755..133S
#
#    The following paper may help to understand general relativistic polarized radiative transfer:
#
#    Shcherbakov R. V., Huang L. 2011, MNRAS, 410, 1052
#              http://adsabs.harvard.edu/abs/2011MNRAS.410.1052S
#
#
#    ASTRORAY is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 3 of the License, or
#    (at your option) any later version.
#
#    ASTRORAY is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ASTRORAY. If not, see <http://www.gnu.org/licenses/>.
#    
#***********************************************************************************/
//Code works fine with double precision accuracy. Experiment with float at your own risk
#define doub double
#include <stdio.h>
#include <iomanip>
#include <string>
#include <math.h>
#include <ctime>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "step.c"
#include "error.c"
#include "cstd.c"
#include "control.c"
#include "evolve.c"
#include "rk2.c"
using namespace std;

//different setups are below
//these files are also platform-specific, as they include directory names in Windows/Linux/UNIX
//#include "win_lin_Jon.c"
#include "win_lin_ADAFs.c"

const doub PI = 4.0*atan(1.0);

const int  ndd=650,           //radial dimension of coordinate/coordinate transformation matrices
	       sflen=14,          //number of frequencies of interest for flux calculations
	       flen=4,            //number of frequencies of interest for images
	       thn=50,            //number of polar angle values to search for 
	       dd=3,              //record size of average temperature & density file
	       wdd=11,            //record size of fluid simulations dump file
	       maxfield=200,      //maximum number of fluid simulations dump files, which can fit in shared memory
	       maxco=3000,        //maximum number of points on a geodesic
	       maxst=12000,       //maximum number of points for radial temperature profile
	       nWlen=120,         //number of frequencies for which propagation coefficients are computed
	       Tlen=100,          //number of temperatures for which propagation coefficients are computed
	       nxy=201,           //actual image resolution in picture plane for imaging (points along a side)
	       snxy=301;          //maximum resolution in picture plane for flux calculations

const doub rgrav=1.33e+12,    //Schwarzschild radius of Sgr A*
	       rrmax=3.4e+5,      //radius in rgrav, where outer temperature and density are defined
	       rhoout=130.,       //outer density for Sgr A*
		   Tout=1.5e+7,       //outer temperature for Sgr A*
		   mp=1.67e-24,       //proton mass
		   me=9.1e-28,        //electron mass
		   cc=3.e+10,         //speed of light
		   kb=1.38e-16,       //Boltzmann constant 
		   ee=4.8e-10,        //electron charge
		   The=me*cc*cc/kb,   //rest mass temperature of electron
		   year=86400.*365.25,//year in seconds
		   Msun=2.00e+33;     //solar mass
const doub r0=20000.;         //maximum radius of each light ray
const doub nWmin=12000.*pow(1.1, -nWlen/2.), nWmax=12000.*pow(1.1, nWlen/2),//minimum and maximum ratios of cyclotron and propagation frequencies, for which propagation effects are non-zero
		   Tmin=0.4, Tmax=0.4*pow(1.05,Tlen), //minimum and maximum ratios of actual and rest mass electron temperatures, for which emissivities are non-zero
		   Tminr=0.4*pow(1.05,-Tlen),         //minimum ratio of actual and rest mass electron temperatures, for which Faraday rotation/conversion are non-zero
		   lnWmin=log(nWmin), lnWmax=log(nWmax), lTminr=log(Tminr), lTmin=log(Tmin), lTmax=log(Tmax);//logarithms

//half-size of the square in a picture plane for each frequency - for flux and image calculations
const doub sftab[sflen][2]={{8.45, 120.}, {14.90, 73.}, {22.50, 63.}, {43.00, 46.}, {87.73, 25.9}, {102., 22.3}, {145., 16.4}, {230.86, 12.2}, {349., 10.3}, {674., 8.8}, {857., 8.6}, {1500., 8.6}, {3000., 8.6}, {5000., 8.6}};

//half-size of the square in a picture plane for each frequency - for image calculations
//const doub freqtab[flen][2]={{43., 45.}, {87., 30.}, {230.86, 12.2}, {690., 7.}};

//polarized spectrum of Sgr A*, each array element is (frequency, Fnu, LP, EVPA, CP)
const doub tofit[sflen][5]={{8.450, 0.683, 0., 0., -0.2500}, {14.90, 0.871, 0., 0., -0.6200}, {22.50, 0.979, 0.1900, 131.0, 0.}, {43.00, 1.135, 0.5500, 94.25, 0.}, {87.73, 1.841, 1.420, -4., 0.}, 
                 {102.0, 1.908, 0., 0., 0.}, {145.0, 2.275, 0., 0., 0.}, {230.9, 2.637, 7.398, 111.5, -1.200}, {349.0, 3.181, 6.499, 146.9, -1.500}, {674.0, 3.286, 0., 0., 0.}, {857.0, 2.867, 0., 0., 0.}, 
                 {1500., 1., 0., 0., 0.}, {3000., 1., 0., 0., 0.}, {5000., 1., 0., 0., 0.}};

//measurement errors of mean fluxes, CP fractions, LP fractions, and EVPAs
const doub dFnu[sflen]={0.031, 0.012, 0.015, 0.026, 0.080, 0.1517, 0.2644, 0.1414, 0.1205, 0.3508, 0.2404, 0., 0., 0.}, //no measurements at highest frequencies
           dCP=0.30, //at 230GHz and 345GHz
		   dLP[3]={0.50, 0.658, 0.605}, //at 87GHz, 230GHz, and 345GHz
		   dEVPA[3]={11.,5.4,2.21};     //at 87GHz, 230GHz, and 345GHz
const bool isLP87=true;//whether to fit for LP fraction at 87GHz. Its observational value is controversial

bool iswrite=true,                                      //whether to write output to a file
	 echeck1=false, echeck2=false, echeck3=false,       //markers for testing (see init.cpp)
	 isBcut=false,                                      //whether to set temperature to zero in certain region close to the BH near polar axis (see evalpointzero.cpp)
	 isBred=false;                                      //whether to reduce temperature in regions of high magnetization (see evalpointzero.cpp)
string fif="";                                          //any modifier for output file name
clock_t start;                                          //timing variable
int fnum,              //fluid simulation dump file number
	cas,               //integer number encoded in LSB_JOBINDEX or PBS_ARRAYID
	sp,                //first command line argument, typically spin
	co,                //second command line argument
	mco,               //third command line argument rounded
	sear,              //fourth command line argument, choice of computation
	ncut,              //the last radial grid point, where the simulation is considered converged
	fdiff=0,           //loading fluid simulation dump files from XXXX-fdiff to XXXX+fdiff to consider simulation evolution as light propagates; fdiff=0 => fast light approximation
	mintim, maxtim,    //physical times of XXXX-fdiff and XXXX+fdiff dump files
	stNt,              //index of electron/ion temperature calculations; is eventually set at 6M radius
	loaded[maxfield];  //numbers of loaded dump files, storing these saves I/O

//Yes, these are global. At least 2 routines use each of those => not always trivial to refactor
doub Bpo,              //third command line argument, often magnetic field strength
     rg,               //horizon radius in units of M
	 fact=1.,          //relative size of integration region, good for tests
	 ans,              //execution time. Too boring to define locally in each place
	 a, asq,           //BH spin and its square
	 th,               //BH spin inclination angle
	 heat,             //electron temperature parameter, determines normalization
	 rhonor,           //density normalization
	 accur,            //relative accuracy of geodesic integration
	 accurr,           //relative accuracy of radiative transfer integration
	 rate,             //accretion rate, typically in g/s
	 minT, maxT,       //minimum and maximum "temperature"=energy density/density
	 ss,               // ss=dr/rg above the BH horizon, where we stop integration of a geodesic
	 fljVc=1., flrQc=1., flrVc=1.,//multiplier to test the behavior of the code for boosted/zeroed V emissivity, Faraday conversion, and Faraday rotation, respectively
	 dphi=0.,          //phi offset to test different phi viewing angles
     TpTe, Te6,        //proton to electron temperature ratio and electron temperature at 6M
	 ts[maxst], te[maxst], tp[maxst], //for computing radial proton and electron temperature profiles
	 rcut,             //radial up to which fluid simulation converged
	 rhopo, rhocon,    //density extension power-law slope and density at rcut
	 Upo, Ucon,        //temperature extension power-law slope and temperature at rcut
	 Bnor,             //magnetic field conversion factor from code units to Gauss
	 rmin,             //inner radius of averaged temperature/density profile (Tsmap*.dat file)
	 thlimit,          //critical parameter for cutting off polar region. Opening angle = arccos(1-thlimit)
	 theta[ndd][thlen],//mapping of code coordinates to physical polar angle theta
	 totin[sflen],  LPo[sflen], CP[sflen], EVPA[sflen], err[sflen],//tofal flux, LP fraction, CP fraction, EVPA, and flux error estimate
	 xtotin[sflen],xLPo[sflen],xCP[sflen],xEVPA[sflen],            //another set of same quantities
     jI[Tlen+1][nWlen+1], jQ[Tlen+1][nWlen+1], jV[Tlen+1][nWlen+1],//emissivities tables
	 rQ[2*Tlen+1], rV[2*Tlen+1],                                   //Faraday conversion and rotation tables
     rtab[2000], Tstab[2000],                                      //for calculating electron and proton temperatures//global since are called in "solvetemperature" routine
     usp[maxfield][phlen][thlen][4], uspKS[maxfield][phlen][thlen];//auxiliary for computing accretion rate in "init" function. Not made local due to stack overflow potential.

float uext[phlen][thlen][5],       //quantities on the spherical fluid simulations convergence boundary
	  dxdxp[ndd][thlen][4][4],     //coordinate transformation smatrix
	  coord[ndd][thlen][2];        //coordinate matrix
typedef float (*uuarr)[phlen][thlen][rlen][wdd]; //type for fluid simulations dump files
              uuarr uu[200];                     //130 dumps fit in 64GB memory for Jon's simulations from 2012

typedef struct {doub lamx[maxco],cooxx[12][maxco];doub llmin,llmax,nu;int indx;} poinx;//geodesic object
     poinx ppy[nthreads]; //define 1 geodesic object per OpenMP thread

//geodesic solver
extern int solvegeodesic(doub t, const doub y[], doub f[], void *params);//line can be commented out
#include "geodes.cpp"

//solving for electron temperature starting from the known outer properties of the system
extern int solvetemperature (doub rz, const doub zz[],doub ff[],void *pas);
#include "solvetemp.cpp"

//initialization
extern int init(int sp, int fmin, int fmax, int sep);
#include "init.cpp"

//polarized radiative transfer
extern int trans (doub llog, const doub yyy[], doub ff[], void *pas);
#include "transnew.cpp"

int main(int argc, char* argv[]) {
	int n1;
    for(n1=0;n1<2*fdiff+1;n1++)
		uu[n1]=(uuarr) new float[phlen][thlen][rlen][wdd];//allocating memory for fluid simulation files
    typedef doub (*ausar)[snxy+1][snxy+1][sflen][5]; 
    ausar ausin = (ausar) new doub[snxy+1][snxy+1][sflen][5];//allocating memory for radiative transfer resuls
    int w,          //thread number for testing
		ittot=0,    //number of steps for radiative transfer solver (sum among all threads)
		niter=0,    //number of iterations for chi^2 minimization
		ix,         //geodesic ID - combined x and y positions
		il,         //auxiliary
		iy,         //geodesic ID along y axis
		kk,         //index of frequency array
		sep=50,     //stepping over fluid simulation snapshot IDs
		kmin=4, 
		kmax=10,    //minimum and maximum frequency array indices of interest
		kstep;      //stepping over frequency array indices
    doub inp[4][3], //variables for chi^2 minimization
		dheat,      //relative variation of heating coefficient
		drho,       //relative vatiation of density
		dtheta;     //absolute variation of (polar) viewing angle of BH spin/disk axis
    start = clock();//timing the entire calculation

    #pragma omp parallel num_threads(nthreads) shared(ittot) private(w)
        {w = omp_get_thread_num();printf("Hello from thread %d\n", w);}//checking that we can create as many threads as we want
    if((descr==NULL) || (argc!=5)){
		printf("Check that 4 numbers are supplied as command line arguments + array ID environment variable is defined\n Exiting\n");
		exit(1);
	}
    sp=atoi(argv[1])-1;  //first command line argument - spin ID
	co=atoi(argv[2]);
	Bpo=atof(argv[3]); 
	mco=floor(Bpo+0.001);//second and third command line arguments are used differently by subroutines
	sear=atoi(argv[4]);  //fourth command line argument - type of computation
//sear=0;//setting a particular type of computation for testing
    switch (sear){
	    case 0: //surf the entire parameter space searching for the best fit to polarized spectrum
	        #include "m_space.cpp"
		    break;
        case 1: //quick computation of spectrum for given sp, heat, rhonor, th
	        #include "m_quick.cpp"
		    break;
        case 2: //convergence studies
	        #include "m_conv.cpp" 
		    break;
        case 3: //surf region close to best fit
	        #include "m_surf.cpp"
		    break;
        case 4: //image of the emitting region
	        #include "m_imag.cpp"
		    break;
        case 5: //search for minimum with a steepest descent method, less reliable than "m_space", but faster
	        #include "m_sear.cpp"
		    break;
        case 6: //averages temperature and density profiles to find Te(Ts) and Tp(Ts) functions
	        #include "m_ts.cpp"
		    break;
	}
    printf ("Iterations per second = % .6f\n", doub(ittot)/ans);
    printf ("Iterations = %d\n ", ittot);
    return(0);
}