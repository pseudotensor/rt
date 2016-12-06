// DEFINE SOME ANSI COLORS FOR MORE READABLE OUTPUT
#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define YELLOW  "\x1b[33m"
#define BLUE    "\x1b[34m"
#define MAGENTA "\x1b[35m"
#define CYAN    "\x1b[36m"
#define RESET   "\x1b[0m"
// UNCOMMENT BELOW TO DEACTIVATE IN CASE IT DOES NOT WORK FOR YOU OR YOU WANT TO REDIRECT OUTPUT TO FILE (YOU WILL GET THINGS LIKE "ESC[33m" IN THERE. SIMPLE TO REMOVE THOUGH...)
// #define RED     ""
// #define GREEN   ""
// #define YELLOW  ""
// #define BLUE    ""
// #define MAGENTA ""
// #define CYAN    ""
// #define RESET   ""

// DEBUG STUFF
int count_nan_gamma = 0;              // ~~~> [setup_avery_toyjet.cpp]
int count_cured_nan_gamma = 0;        // ~~~> [setup_avery_toyjet.cpp]
doub smallest_radius_where_nan = 1e9; // ~~~> [setup_avery_toyjet.cpp]
doub largest_radius_where_nan  = 0.;  // ~~~> [setup_avery_toyjet.cpp]

// GEODESIC DIAGNOSTIC
// #DEFINE GEODESIC_DIAGNOSTIC true ?
bool GEODESIC_DIAGNOSTIC=true; //true; // false;

int geodesic_output_x=0; // choose pixel to focus on
int geodesic_output_y=0; // choose pixel to focus on
int geodesic_output_every_x = 5; // output geodesic information for every nth geodesic along x-dir in image plane
int geodesic_output_every_y = 5; // output geodesic information for every nth geodesic along y-dir in image plane

bool TEMPERATURE_DIAGNOSTIC=false; // output info on t_p,t_e,u,rho,r,th,ph
string TEMPERATURE_PRESCRIPTION="sharma+isoth"; // sharma , sharma+isoth , constant_tetp_fraction
// string TEMPERATURE_PRESCRIPTION="constant_tetp_fraction"; // "sharma+isoth"; // sharma , sharma+isoth , constant_tetp_fraction
// string TEMPERATURE_PRESCRIPTION="sharma"; // "sharma+isoth"; // sharma , sharma+isoth , constant_tetp_fraction
const string INTERPOLATION_METHOD="ASTRORAYv1.0"; // "SUM"; // can be "SUM","ASTRORAYv1.0",...
const string DRIVE_FIT_TO="ASTRORAYv1.0"; // can be "DEXTER","CK","ASTRORAYv1.0","ASTRORAYv1.0_unpolarized",...

bool BREMSSTRAHLUNG=false; // INCLUDE ee & ei BREMSSTRAHLUNG EMISSIVITY ? FIXME:VALUES SEEM WRONG NEED2DEBUG 

const int nr_of_imagediags=5; // 4:IQUV + images of other quantities  see [intensity.cpp] and [solvetrans.cpp]

const doub PI = 4.0*atan(1.0);

const bool avoid_pole=true; // grep for critan ~> [evalpointzero.cpp]
// const bool avoid_pole=false; // grep for critan ~> [evalpointzero.cpp] // gives trouble

//Sgr A see math/checks_GRMHD_code.nb->Intensity functions->Check definitions of distance to Sgr A*, percentages, position angles
doub Jy2cgs        = 1e23; // 1Jy=10^23 erg/(s cm^2 Hz) [cgs]
#if SOURCE==SgrA
const doub rgrav=1.33e+12;    //Schwarzschild radius of Sgr A* //RG: in cm corresponds to M_BH~4.4e6 Msun RG:TODO RENAME TO rs!
// const doub rgrav=2.66e+12;    // Gravitational radius of Sgr A* //RG: in cm corresponds to M_BH~4.4e6 Msun RG:TODO RENAME TO rs!
doub ang_size_norm = 66.4648/Jy2cgs; // 10^23 (rg^2/D^2) // D=8.4kpc // M=4.5e6Msun 
// doub ang_size_norm = 196.5443582815904/Jy2cgs; // 10^23 (rg^2/D^2) // D=8.4kpc=3e22cm // M=4.5e6Msun 
#elif SOURCE==M87
// const doub rgrav=1.33e+15;    //Schwarzschild radius of M87 //RG: in cm corresponds to M_BH~4.4e9 Msun RG:TODO RENAME TO rs!
const doub rgrav=1.099e+15;    //Schwarzschild radius of M87 //RG: in cm corresponds to M_BH~3.4e9 Msun RG:TODO RENAME TO rs! // Broderick & Loeb 2009
doub ang_size_norm = pow(rgrav/1.33e12,2)*pow(8.4/16.7e3,2)*66.4648/Jy2cgs; // 10^23 (rg^2/D^2) // D=?Mpc // M=?Msun 
#endif 
// MODELS
// const char avery_toy_jet[64]="yes"; // global flag  to turn on/off Avery's toyjet + RIAF model 
const char avery_toy_jet[64]="no"; // global flag  to turn on/off Avery's toyjet + RIAF model 

// bool use_radial_extension=false; //RG: SET rho=0 (~> no emission/absorption/FR/FC) outside rcut
bool use_radial_extension=true; //RG: USE radial extension outside rcut
const string RADIAL_EXT_FILE="dxdxp.dat"; // "usgdump2d";

const string image_diagnostic="lambda"; // "melrose","column density","lambda" see [transnew.cpp]

const int sflen=14,          //number of frequencies of interest for flux calculations
  flen=4,            //number of frequencies of interest for images
  thn=50,            //number of polar angle values to search for 
  dd=3,              //record size of average temperature & density file

// THICKDISK7,THICKDISKHR3?,DIPOLE3DFIDUCIALA,QUADRUPOLE
#if MODEL==THICKDISK7 || MODEL==DIPOLE || MODEL==QUADRUPOLE
  wdd=11,            //record size of fluid simulations dump file
// a=0 MAD rtf2_15r35_a0.0_0_0_0 , thinnermad*
#elif MODEL==THINNERMAD
  wdd=11+3,          //record size of fluid simulations dump file
#endif
  maxfield=200,      //maximum number of fluid simulations dump files, which can fit in shared memory
//maxco=3000,        //maximum number of points on a geodesic
// on bh01 maxco=80000 is ok, but maxco=90000 the code just silently quits ~> exceed heap size limit??
//  maxco=50000,        //maximum number of points on a geodesic
  maxco=80000,        //maximum number of points on a geodesic
  maxst=40000,       //maximum number of points for radial temperature profile
  nWlen=120,nWlen_nth=120,         // number of frequency bins for lookup tables of propagation coefficients // nWlen=60 
  Tlen=100,Tlen_nth=160/*160*/,          // number of temperature bins for lookup tables of propagation coefficients 
  nxy=100, // 1, // 299, //           //actual image resolution in picture plane for imaging (points along a side)
  snxy=nxy;          //maximum resolution in picture plane for flux calculations

// doub rrmax=3.4e5;
// doub rrmax=95.477943855303153;

const doub 
           // DEFAULT
           rrmax=3.4e+5,      //radius in rgrav, where outer temperature and density are defined
           // WIP
// rrmax=7.11445236e+03,      //radius in rgrav, where outer temperature and density are defined
// rrmax=668.39419037309119,      //radius in rgrav, where outer temperature and density are defined
	       rhoout=130.,       //outer density for Sgr A*
		   Tout=1.5e+7,       //outer temperature for Sgr A*
		   mp=1.67e-24,       //proton mass [g]
		   me=9.1e-28,        //electron mass [g]
		   cc=3.e+10,         //speed of light
		   kb=1.38e-16,       //Boltzmann constant 
		   ee=4.8e-10,        //electron charge
           h_planck=6.6260755e-27, // Planck's constant h (not hbar) units: [erg s]
		   The=me*cc*cc/kb,   //rest mass temperature of electron
		   year=86400.*365.25,//year in seconds
           Msun=2.00e+33,     //solar mass [g]

// r0=20000.;         //maximum radius of each light ray: DEFAULT ASTRORAY v1.0
// r0=500.;         //maximum radius of each light ray // RG: coordinate distance (in rgrav units) of image plane to BH (Horizon)?
// r0=1000.;         //maximum radius of each light ray 
r0 = 668.39419037309119; // for grtrans comparison

// Temperature sampling & range for propagation effects for THERMAL
const doub nWmin=12000.*pow(1.1, -nWlen/2.), nWmax=12000.*pow(1.1, nWlen/2),//minimum and maximum ratios of cyclotron and propagation frequencies, for which propagation effects are non-zero
		   Tmin=0.4, Tmax=0.4*pow(1.05,Tlen), //minimum and maximum ratios of actual and rest mass electron temperatures, for which emissivities are non-zero
		   Tminr=0.4*pow(1.05,-Tlen),         //minimum ratio of actual and rest mass electron temperatures, for which Faraday rotation/conversion are non-zero
		   lnWmin=log(nWmin), lnWmax=log(nWmax), lTminr=log(Tminr), lTmin=log(Tmin), lTmax=log(Tmax);//logarithms

// Temperature sampling & range for propagation effects for NON-THERMAL
const doub logspacing_Wmin_nth=1.25, logspacing_Wmax_nth=1.25;
const doub logspacing_Tmin_nth=1.05, logspacing_Tmax_nth=1.05;
const doub nWmin_nth=12000.*pow(logspacing_Wmin_nth, -nWlen_nth/2.), nWmax_nth=12000.*pow(logspacing_Wmax_nth, nWlen_nth/2),//minimum and maximum ratios of cyclotron and propagation frequencies, for which propagation effects are non-zero
		   Tmin_nth=0.4, Tmax_nth=0.4*pow(logspacing_Tmin_nth,Tlen_nth), //minimum and maximum ratios of actual and rest mass electron temperatures, for which emissivities are non-zero
		   Tminr_nth=0.4*pow(logspacing_Tmax_nth,-Tlen_nth),         //minimum ratio of actual and rest mass electron temperatures, for which Faraday rotation/conversion are non-zero
		   lnWmin_nth=log(nWmin_nth), lnWmax_nth=log(nWmax_nth), lTminr_nth=log(Tminr_nth), lTmin_nth=log(Tmin_nth), lTmax_nth=log(Tmax_nth);//logarithms


//half-size of the square in a picture plane for each frequency - for flux and image calculations
// RG: frequencies are: 8.45, 14.90, 22.50, 43.00, 87.73, 102., 145., 230.86, 349., 674., 857., 1500., 3000., 5000. GHz?
// RG: {{frequency1, half-screen-size@frequency1?}, {frequency2, half-screen-size@frequency2?}, ...}

//DEFAULT (SGR A*)
// const doub sftab[sflen][2]={{8.45, 120.}, {14.90, 73.}, {22.50, 63.}, {43.00, 46.}, {87.73, 25.9}, {102., 22.3}, {145., 16.4}, {230.86, 12.2}, {349., 10.3}, {674., 8.8}, {857., 8.6}, {1500., 8.6}, {3000., 8.6}, {5000., 8.6}};

// const doub bandwidth[sflen][2]={{8.45, 120.}, {14.90, 73.}, {22.50, 63.}, {43.00, 46.}, {87.73, 25.9}, {102., 22.3}, {145., 16.4}, {230.86, 15.0}, {349., 10.3}, {674., 8.8}, {857., 8.6}, {1500., 8.6}, {3000., 8.6}, {5000., 8.6}};
doub bandwidth=0.; //15.; // 0.; // choose non-zero to output sideband
const doub sftab[sflen][2]={{8.45, 120.}, {14.90, 73.}, {22.50, 63.}, {43.00, 46.}, {87.73, 25.9}, {102., 22.3}, {145., 16.4}, {230.86 + bandwidth, 15.0}, {349., 10.3}, {674., 8.8}, {857., 8.6}, {1500., 8.6}, {3000., 8.6}, {5000., 8.6}};

// DEXTER 40Mx40M @ 230GHz
// const doub sftab[sflen][2]={{8.45, 120.}, {14.90, 73.}, {22.50, 63.}, {43.00, 46.}, {87.73, 25.9}, {102., 22.3}, {145., 16.4}, {230.86, 10.0}, {349., 10.3}, {674., 8.8}, {857., 8.6}, {1500., 8.6}, {3000., 8.6}, {5000., 8.6}};
//M87 WIP: STARTED MODIFYING 230Ghz entry... but bizarre code behavior...
//const doub sftab[sflen][2]={{8.45, 120.}, {14.90, 73.}, {22.50, 63.}, {43.00, 46.}, {87.73, 25.9}, {102., 22.3}, {145., 16.4}, {230.86, 5.0}, {349., 10.3}, {674., 8.8}, {857., 8.6}, {1500., 8.6}, {3000., 8.6}, {5000., 8.6}};

// polarized spectrum of Sgr A*, each array element is 
// {frequency, Fnu, LP, EVPA, CP}
// Fnu:  flux at frequency nu
// LP:   linear polarization fraction (image-averaged/zero-baseline) at frequency nu
// EVPA: Electric vector position angle (image-averaged/zero-baseline)
// CP:   circular polarization fraction (image-averaged/zero-baseline) at frequency nu 
const doub tofit[sflen][5]={
  {8.450, 0.683, 0., 0., -0.2500}, 
  {14.90, 0.871, 0., 0., -0.6200}, 
  {22.50, 0.979, 0.1900, 131.0, 0.}, 
  {43.00, 1.135, 0.5500, 94.25, 0.}, 
  {87.73, 1.841, 1.420, -4., 0.}, 
  {102.0, 1.908, 0., 0., 0.}, 
  {145.0, 2.275, 0., 0., 0.}, 
  {230.9, 2.637, 7.398, 111.5, -1.200}, 
  {349.0, 3.181, 6.499, 146.9, -1.500}, 
  {674.0, 3.286, 0., 0., 0.},
  {857.0, 2.867, 0., 0., 0.},
  {1500., 1., 0., 0., 0.},{3000., 1., 0., 0., 0.},{5000., 1., 0., 0., 0.}};

// measurement errors of mean fluxes, CP fractions, LP fractions, and EVPAs (no measurements at highest frequencies)
const doub dFnu[sflen]={0.031, 0.012, 0.015, 0.026, 0.080, 0.1517, 0.2644, 0.1414, 0.1205, 0.3508, 0.2404, 0., 0., 0.}, 
           dCP=0.30,                    // at 230GHz and 345GHz
           dLP[3]={0.50, 0.658, 0.605}, // at 87GHz, 230GHz, and 345GHz
           dEVPA[3]={11.,5.4,2.21};     // at 87GHz, 230GHz, and 345GHz
const bool trustLP87=true;//whether to fit for LP fraction at 87GHz. Its observational value is controversial
doub dof=7.;                            //degrees of freedom

doub heat_min=0.3,                // range of heat constants probed in [m_space.cpp]
  heat_max=0.75;

bool nth=false,                                          // include non-thermal electrons?
//bool nth=true,                                          // include non-thermal electrons?

     iswrite=true,                                      //whether to write output to a file
	 echeck1=false, echeck2=false, echeck3=false,       //markers for testing (see init.cpp)
	 isBcut=false,                                      //whether to set temperature to zero in certain region close to the BH near polar axis (see evalpointzero.cpp)
	 isBred=false;                                      //whether to reduce temperature in regions of high magnetization (see evalpointzero.cpp)
doub include_jet = 1.;
doub Te_jet_par = 10.; // isothermal jet temperature in [me*cc*cc/kb]
doub magn_cap=10.; //used to reduce temperature/rho/emission/absorption in
doub magn_cap_rho=4.; //used to reduce temperature/rho/emission/absorption in regions of high magnetization (see evalpointzero.cpp)
doub magn_floor=0.; //used to lighten-up the jet and turn off the disk (see [evalpointzero.cpp])
doub trace_theta_slice_width=10.*PI/100.; // used to trace geodesics arising from a thin slice in the theta direction. Reduce temperature/rho/emission/absorption in all other regions [see evalpointzero.cpp]
doub trace_theta_slice_angle; // used to trace geodesics arising from a thin slice in the theta direction. Reduce temperature/rho/emission/absorption in all other regions [see evalpointzero.cpp]
doub trace_r_slice_width=2.; // used to trace geodesics arising from a thin slice in the theta direction. Reduce temperature/rho/emission/absorption in all other regions [see evalpointzero.cpp]
doub trace_r_slice; // used to trace geodesics arising from a thin slice in the r direction. Reduce temperature/rho/emission/absorption in all other regions [see evalpointzero.cpp]
string fif="";                                          //any modifier for output file name
clock_t start;                                          //timing variable
int fnum,              //fluid simulation dump file number
	cas,               //integer number encoded in LSB_JOBINDEX or PBS_ARRAYID

/*********************************/
/********* CMD-LINE-ARGS *********/
	sp,                //first command line argument, typically spin
	co,                //second command line argument
	mco,               //third command line argument rounded
	sear,              //fourth command line argument, choice of computation
/********* CMD-LINE-ARGS *********/
/*********************************/

	ncut,              //the last radial grid point, where the simulation is considered to be in the final/relaxed/asymptotic state
	fdiff=200,           //loading fluid simulation dump files from XXXX-fdiff to XXXX+fdiff to consider simulation evolution as light propagates; fdiff=0 => fast light approximation
	mintim, maxtim,    //physical times of XXXX-fdiff and XXXX+fdiff dump files
	stNt,              //index of electron/ion temperature calculations; is eventually set at 6M radius 
	loaded[maxfield];  //numbers of loaded dump files, storing these saves I/O


//Yes, these are (RG: non-constant) global. At least 2 routines use each of those => not always trivial to refactor
doub Bpo,              //third command line argument, often magnetic field strength
     rg,               //horizon radius in units of M
//fact=1.,          //relative size of integration region, good for tests
	 fact=1.,          //relative size of integration region, good for tests
	 ans,              //execution time. 
	 a, asq,           //BH spin and its square
	 th,               //BH spin inclination angle // RG: Do we mean latitude/inclination angle w.r.t. BH spin?!
	 heat,             //electron temperature parameter, determines normalization // same as "C" in [Shcherbakov,Penna,McKinney 2012] eq ??
	 rhonor,           //density normalization/unit/scale
	 accur,            //relative accuracy of geodesic integration
	 accurr,           //relative accuracy of radiative transfer integration
	 rate,             //accretion rate, typically in g/s
	 minT, maxT,       //minimum and maximum "temperature"=energy density/density
	 ss,               // ss=dr/rg above the BH horizon, where we stop integration of a geodesic
	 flrQc=1., flrVc=1.,//multiplier to test the behavior of the code for boosted/zeroed Faraday conversion, and Faraday rotation, respectively
	 fljIc=1., fljQc=1., fljVc=1.,//multiplier to test the behavior of the code for boosted/zeroed I,Q,V emissivities
	 dphi=0.,          //phi offset to test different phi viewing angles
  Te_jet, TpTe_jet, TpTe, Te6,        //proton to electron temperature ratio and electron temperature at 6M
	 ts[maxst], te[maxst], tp[maxst], //for computing radial proton and electron temperature profiles
	 rcut,             //radius up to which fluid simulation converged // RG: converged->relaxed to steady state?
	 rhopo, rhocon,    //density extension power-law slope and density at rcut
	 Upo, Ucon,        //temperature extension power-law slope and temperature at rcut
	 Bnor,             //magnetic field conversion factor from code units to Gauss
	 rmin,             //inner radius of averaged temperature/density profile (Tsmap*.dat file)
	 thlimit,          //critical parameter for cutting off polar region. Opening angle = arccos(1-thlimit)
	 theta[ndd][thlen],//mapping of code coordinates to physical polar angle theta
	 totin[sflen],  LPo[sflen], CP[sflen], EVPA[sflen], err[sflen],//total flux, LP fraction, CP fraction, EVPA, and flux error estimate
	 xtotin[sflen],xLPo[sflen],xCP[sflen],xEVPA[sflen],            //another set of same quantities // RG: but not copy of the former... so explain the difference! Dont use x!
     jI[Tlen+1][nWlen+1], jQ[Tlen+1][nWlen+1], jV[Tlen+1][nWlen+1],//emissivities tables THERMAL
	 rQ[2*Tlen+1], rV[2*Tlen+1],                                   //Faraday conversion and rotation tables THERMAL
     jI_nth[Tlen_nth+1][nWlen_nth+1], jQ_nth[Tlen_nth+1][nWlen_nth+1], jV_nth[Tlen_nth+1][nWlen_nth+1],//emissivities tables NON-THERMAL
     aI_nth[Tlen_nth+1][nWlen_nth+1], aQ_nth[Tlen_nth+1][nWlen_nth+1], aV_nth[Tlen_nth+1][nWlen_nth+1],//absorptivity tables NON-THERMAL
	 rQ_nth[2*Tlen_nth+1], rV_nth[2*Tlen_nth+1],                           //Faraday conversion and rotation tables NON-THERMAL
     rtab[2000], T_sim_tab[2000],                                      //for calculating electron and proton temperatures//global since are called in "solvetemperature" routine
     usp[maxfield][phlen][thlen][4], uspKS[maxfield][phlen][thlen];//auxiliary for computing accretion rate in "init" function. Not made local due to stack overflow potential.

float uext[phlen][thlen][5],       //quantities on the spherical fluid simulations convergence boundary
	  dxdxp[ndd][thlen][4][4],     //coordinate transformation smatrix
	  coord[ndd][thlen][2];        //coordinate matrix
typedef float (*uuarr)[phlen][thlen][rlen][wdd]; //type for fluid simulations dump files

      uuarr uu[200];                     //130 dumps fit in 64GB memory for Jon's simulations from 2012
