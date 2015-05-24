/***********************************************************************************
#    Copyright 2014 Roman Shcherbakov
#                   
#              ASTRORAY version 1.0   (released June 25, 2014)
#
#    ASTRORAY v1.0 is a program that performs general relativistic polarized radiative transfer 
#    of synchrotron radiation near black holes. The code employs a ray tracing technique 
#    and can easily handle large optical depth (~10,000). 
#    ASTRORAY produces images as well as spectra.
#    
#    The latest version of ASTRORAY code can be retrieved from 
#    our distribution website: http://astroman.org/code/ASTRORAY/.
#
#    This version of ASTRORAY is configured to use input files from a 
#    3D version of HARM code maintained by Jonathan McKinney (UMD).
#    The code assumes that the source is a plasma near a
#    black hole described in Kerr-Schild coordinates.
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
#define doub double // PRECISION: double: ok. float: Use at your own risk!
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

#include <unistd.h>       // -> sleep
#include <sys/resource.h> // -> getrusage -> PROFILING

using namespace std;

// different setups are below
// these files are also platform-specific, as they include directory names in Windows/Linux/UNIX
#include "win_lin_Jon.c" // #include "win_lin_ADAFs.c"

#include "profiling.cpp"
#include "global_variables.cpp" // contains most global variables

typedef struct {doub lamx[maxco],cooxx[12][maxco];doub llmin,llmax,nu;int indx;} poinx;//geodesic object
     poinx ppy[nthreads]; //define 1 geodesic object per OpenMP thread

//geodesic solver
extern int solvegeodesic(doub t, const doub y[], doub f[], void *params);//line can be commented out
#include "geodes.cpp"

//solving for electron temperature starting from the known outer properties of the system
extern int solvetemperature (doub rz, const doub zz[],doub ff[],void *pas);
#include "solvetemp.cpp"

//RG:CLEANUP!
//initialization
//RG:FIXME "UNDEFINED REFERENCE"
//extern int setup_averys_toyjet(int i,int j,int k, float uu[phlen][thlen][rlen][wdd]);//, float Br,float Btheta,float Bphi,float ur,float utheta,float uphi,float rhoL, doub restL);
//extern int setup_averys_toyjet(int i,int j,int k);//, float Br,float Btheta,float Bphi,float ur,float utheta,float uphi,float rhoL, doub restL);
//extern int setup_averys_toyjet(doub coord, float Br,float Btheta,float Bphi,float ur,float utheta,float uphi,float rhoL, doub restL);
//extern int setup_averys_toyjet();//, float restL);
extern int init(int sp, int fmin, int fmax, int sep); //RG:CLEANUP/INVESTIGATE sep not used
#include "init.cpp"

//polarized radiative transfer
extern int trans (doub llog, const doub yyy[], doub ff[], void *pas);
//extern int trans (doub llog, const doub yyy[], doub ff[], void *pas, int stNx); //RG:WIP
#include "transnew.cpp"



/**********************************
 ********* MAIN ROUTINE ***********
 **********************************/

int main(int argc, char* argv[]) {

	int n1;
    struct rusage usage; // FOR PROFILING PURPOSES
    
    // RG:FIXME BAD PLACE TO DO THAT HERE. fdiff can change (e.g. in cases)
    // RG: -> should make fdiff const global!
    for(n1=0;n1<2*fdiff+1;n1++)
      uu[n1]=(uuarr) new float[phlen][thlen][rlen][wdd]; // allocate memory for fluid data 

    typedef doub (*ausar)[snxy+1][snxy+1][sflen][5]; 
    ausar ausin = (ausar) new doub[snxy+1][snxy+1][sflen][5]; // allocate memory for radiative transfer results
    int w,          //thread number (for testing)
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
    start = clock();//timing entire computation //RG: ADD MORE BENCHMARKING/PROFILING see [profiling.cpp]

    #pragma omp parallel num_threads(nthreads) shared(ittot) private(w)
    { // checking that we can create as many threads as we want
      w = omp_get_thread_num();
      // printf(YELLOW"[ASTRORAY_main.cpp]:"RESET" Hello from thread %d\n", w);
    }

    if((descr==NULL) || (argc!=5)){
		printf(YELLOW"[ASTRORAY_main.cpp]:"RESET"  Check that 4 numbers are supplied as command line arguments + array ID environment variable is defined\n Exiting\n");
		exit(1);
	}


/*********************************/
/********* CMD-LINE-ARGS *********/
    sp=atoi(argv[1])-1;   //first command line argument - index for spin ID
	co=atoi(argv[2]);     // argv[4]=6:fmin [m_ts.cpp]
	Bpo=atof(argv[3]);    // 
	//mco=floor(Bpo+0.001); // argv[4]=6:fmax [m_ts.cpp] (second and third command line arguments are used differently by subroutines)
	mco=atoi(argv[3]); // argv[4]=6:fmax [m_ts.cpp] (second and third command line arguments are used differently by subroutines)
	sear=atoi(argv[4]);  //fourth command line argument - type of computation
/********* CMD-LINE-ARGS *********/
/*********************************/


    /**** CHOOSE MODE OF OPERATION ****/
    switch (sear){
	    case 0: //Scan full parameter space for best fit to polarized spectrum
            cout << YELLOW"[ASTRORAY_main.cpp]:"RESET" ENTERING m_space..."<<endl;
	        #include "m_space.cpp"
		    break;
        case 1: //quick computation of spectrum for given sp, heat, rhonor, th
            cout << YELLOW"[ASTRORAY_main.cpp]:"RESET" ENTERING m_quick..."<<endl;
	        #include "m_quick.cpp"
		    break;
        case 2: //convergence studies //RG:sensitivity to various parameters
            cout << YELLOW"[ASTRORAY_main.cpp]:"RESET" ENTERING m_conv..."<<endl;
	        #include "m_conv.cpp" 
		    break;
        case 3: //surf region close to best fit //RG:surf="search in more detail"
            cout << YELLOW"[ASTRORAY_main.cpp]:"RESET" ENTERING m_surf..."<<endl;
	        #include "m_surf.cpp"
		    break;
        case 4: //image of the emitting region
            cout << YELLOW"[ASTRORAY_main.cpp]:"RESET" ENTERING m_imag..."<<endl;
	        #include "m_imag.cpp"
		    break;
        case 5: //search for minimum with a steepest descent method, less reliable than "m_space", but faster
            cout << YELLOW"[ASTRORAY_main.cpp]:"RESET" ENTERING m_sear..."<<endl;
	        #include "m_sear.cpp"
		    break;
        case 6: //averages temperature and density profiles to find Te(Ts) and Tp(Ts) functions
            cout << YELLOW"[ASTRORAY_main.cpp]:"RESET" ENTERING m_ts..."<<endl;
	        #include "m_ts.cpp"
		    break;
	}

    printf(YELLOW"[ASTRORAY_main.cpp]"BLUE"==============PROFILING INFO==============\n"RESET);
    // Use gprof to obtain time taken by each function. Use these flags for compilation -pg -fprofile-arcs -ftest-coverage.

    // RG: BELOW PRINTS 0.01 instead of 1
    // time_t t_test=clock();
    // sleep(1);
    // printf(YELLOW"[ASTRORAY_main.cpp]:"RESET" TEST TIMER: %f",(clock()-t_test)/(doub)CLOCKS_PER_SEC);
    doub t_fudge=0.01; // test with sleep(1)
    doub t_total = (clock()-start)/(doub)CLOCKS_PER_SEC;
    printf(YELLOW"[ASTRORAY_main.cpp]"RESET" TOTAL RUNTIME             = %.1f secs (%.1f secs across all threads)\n", t_total* t_fudge, t_total *t_fudge * nthreads);
    //printf ("===================RT=====================\n");
    //printf ("RUNTIME RT SOLVER         = % .1f secs (%.0f%% %.1f across all threads)\n", ans*t_fudge, ans/t_total*nthreads, ans*t_fudge);
    printf(YELLOW"[ASTRORAY_main.cpp]"RESET" Total iterations          = %d\n", ittot);
    printf(YELLOW"[ASTRORAY_main.cpp]"RESET" Computational performance = %.0f iterations/sec\n", doub(ittot)/t_total/t_fudge);
    printf(YELLOW"[ASTRORAY_main.cpp]"RESET" RUNTIME GEODESICS         = %.1f secs (%.0f%%)\n", t_geodesics *t_fudge, t_geodesics/t_total/nthreads*100);
    printf(YELLOW"[ASTRORAY_main.cpp]"RESET" RUNTIME SOLVETRANS        = %.1f secs (%.0f%%)\n", t_solvetrans*t_fudge, t_solvetrans/t_total/nthreads*100);
    printf(YELLOW"[ASTRORAY_main.cpp]"BLUE"============PROFILING INFO END============\n"RESET);



    cout << GREEN"====\nDONE\n===="RESET << endl;

    return(0);
}
