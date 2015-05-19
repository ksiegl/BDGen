#define PI 3.141592653589793238
#define ALPHA 0.000000000000076375		/* value from NIM A 404 (1998) 305-310 */
#define mexp 0.51099906         /* mexp and merr in MeV */
#define merr 0.00000015

#define GAMMANUM 512
#define GTNUM 1024
#define LVLNUM 512

#define beta_states 512
#define neutron_states 512
#define charge_states 16
#define fermi_points 49  // This is constant for all isotopes
#define Egammaexp2 0.000000	// gamma energy and error in MeV 
#define Egammaerr2 0.000011

#define Eshakeexp  0.000000	/* guess 10 eV ion energy */
#define uexp 931.494028          /* conversion between amu and MeV */
#define uerr   0.000023

#define M_n_exp      1008664.91574    /* neutron rest mass in 1e-6 amu from AME2003 */
#define M_n_err            0.00056

//#define Sn_exp              0.0
//#define Sn_exp             7.675     /* neutron separation energy of Te-137 from AME2012*/
#define Sn_err             0.005     /* unc from AME2012 */

//#define E_n_exp            0.2        /* neutron energy in MeV */
#define E_n_err            0.001        /* I don't know this uncertainty */

#define E_two_alpha_exp    0.092        /* energy of 8Be above 2 alpha */

#define ftexp 4106.4            /* 21Na ft-value (seconds) */
#define fterr   11.6
#define SIDEBR		0.00019		/* the 0+ --> 0+ branching ratio */
#define HARDBR	  	0.02518		/* contribution from hard brem. -- must be calculated -- 0.02518 */
#define SIDEHARDBR	0.02243		/* contribution to side branch from hard brem. 	      -- 0.02243 */
#define CS	0.001		/* see F. Gluck/Computer Phys. Comm. 101 (1997) 223-231. */
#define BREMMAX	0.065		/* max of the brem phase space -- must be figured out -- 0.065 */

#define trap_x     0.0      /* should be 0 */
#define trap_y    0.0      /* center is 0   */
#define trap_z    0.0      /* center is 0   */
#define CENTER_X   0.0    /* center of chamber is 0 */
#define CENTER_Y   0.0	/* center of chamber is 0 */
#define CENTER_Z  0.0	/* center of chamber is 0 */
extern double TRAP_FWHM;
#define TRAP_FWHM_XSCALE 1.0
#define TRAP_FWHM_YSCALE 1.0
#define TRAP_FWHM_ZSCALE 1.0

//#define CHARGE 2
#define RF_FREQ 0.310  /* rf frequency in MHz */

// Commenting out old values of BD_X and BD_Y
//#define BD_X -36.0      /* Beta Detector distance from center of trap, used to be 55.16 in POP */
//#define BD_Y 75.0      /* Beta Detector distance from center of trap, used to be 61.3 in POP */
#define BD_X 0.0	  	/* Beta Detector distance from center of trap */
#define BD_Y 0.0 		
#define BD_Z 85.3      /* Beta Detector distance from center of trap, used to be 60 in POP */
#define BD_SIDE 50.8 /* length of side of Si (each is 50.5 with 16 strips) */ 
#define NUM_STRIPS 16.0 /* number of strips on the strip detector */     
#define BD_R_MAX  3     /* next integer above BD_RADIUS+2 */
#define BD_RADIUS 41.148  /* radius of Beta Detector in mm */

#define GE_DIST	82.55	/* distance of HPGe from trap center */
#define GE_R    38.05           /* HPGe crystal radius */
#define GAMMA_ATTEN 0.0	/* variation in gamma-ray detection efficiency between edge and center */
#define NUM_GE	4	/* number of HPGe detectors */

#define NUMCHANNELS	200000
#define SIGMA	0.00972
#define COMPTON	0.667

#define ACCLENGTH 254	/* number of divisions of acceptance.txt file */
#define BEFILEMAX 1010	/* number of bins of the Be window energy loss file -- 1010 */
#define BEFILESCALE 0.0025	/* scale for Be window loss file (in MeV) */
#define E_FIELD_LOSS 0.005	/* positron energy lost due to electric field (in MeV) */

#define DIVISIONS 10000
#define VSACCURACY 100	/* number of terms in sum of Spence function */

#define ALPHA_ENERGY_FWHM 0.0000000000000000078 /* FWHM of Si detector for alphas in units of me 0.078 => 40 keV */
#define BETA_ENERGY_FWHM 0.0000000000000015 /* fractional uncertainty at energy of me (511 keV) and scales as sqrt(E_beta) */
#define BETA_THRESHOLD 0.0 /* threshold for beta detection in units of me, originally 0.2446=125keV or .0587=30keV */
#define GE_FWHM 	0.00000000000000035	/* FWHM of Ge resolution in units of me 0.0035 => 1.8 keV */
#define SI_FWHM		0.00972			/* testing with Si resolution of FWHM of 500 keV */

#define NFILETYPES 1
#define NUMSPEC 1	/* number of different files... _X, _Xhard, _Xside, _Xsidehard */
//#define a      -0.333	/* value the Standard Model predicts for main branch, was 0.333 */
#define MF	0.000	/* Holstein's "a", or Fermi strength */
#define MGT	1.000	/* Holstein's "c", or Gamow-Teller strength */
#define bWM	82.63	/* Holstein's "b", or weak magnetism */
#define dSCC	0.000	/* Holstein's "d", or induced tensor, second-class current */
#define bFierz	0.000	/* Fierz interference term */
#define a_side 	       -0.333   /* value the Standard Model predicts for side branch, including corrections */
#define A_side         -0.600   /* beta-asymmetry coefficient */
#define B_side       	0.600   /* neutrino-asymmetry coefficient */
#define c_align_side   -0.200   /* alignment term */
#define MF_side		0.000
#define MGT_side	0.391
#define bWM_side	21.8
#define dSCC_side	0.000
#define bFierz_side	0.000
#define Dtrv_side	0.000
#define Z		8.0
#define NUCRAD		0.01
#define SpeedofLight    2.9979e11   /* speed of light in mm/sec */

#define beta_neutron_p    0.0   /* just a guess for this correlation term */

/* only needed if a polarized sample were to be used */
#define A		0.000	/* beta-asymmetry term */
#define B		0.000	/* neutrino-asymmetry term */
#define c_align		0.000	/* alignment term */
#define Dtrv		0.000	/* time-reversal-violating D coefficient */
static double polarization[3]={0.000,0.000,0.000};	/* nuclear polarization in x,y,z */
static double alignment[3]={0.000,0.000,0.000};		/* nuclear alignment in x,y,z */
#define ALIGN		0.000
/* end polarized stuff */

static char GAMMA_SIMION_FILE[40]= "gamma_simion_input.ion";
static char GAMMA_GEANT_FILE[40]="gamma_geant_input.txt";
static char GAMMA_AUXFILE[40]="gamma_auxfile.txt";
static char NEUTRON_SIMION_FILE[40]= "neutron_simion_input.ion";
static char NEUTRON_GEANT_FILE[40]= "neutron_geant_input.txt";
static char NEUTRON_AUXFILE[40]= "neutron_auxfile.txt";

void message1(int file);

#define NUM_BE_BINS 20
#define BEHIST_START 1.00
#define BEHIST_END 5.00
#define NUM_GE_BINS 100
#define GEHIST_START 4.5204
#define GEHIST_END	4.5304

/* these are the files used in BDSimFile.c */
#define DATAFILE_1 "BETA1.mac"
#define DATAFILE_A "BETA_1.mac"
#define OUTPUT_1 "TOFData_1"
#define OUTPUT_A "TOFData_a"
#define HIST_1 "TOFHist_1"
#define HIST_A "TOFHist_a"
#define AVE_1 "E_betaHist_1"
#define AVE_A "E_betaHist_a"

/*	The energylevel struct stores the energy, number of gamma decays, and total transition normalization
*	for each level listed in the RIPL-3 file.  firstgamma is the index of the first of the listed transitions
*	for the level.  All of the transitions for a level are indexed firstgamma,...,firstgamma+numofgamma-1
*/
typedef struct energylevel
{
	double energy;
	int numofgamma;
	int firstgamma;
	double totprob;
}ELevel;

/*	Each gammatransition contains the index of its initial state, ilevel, its final state, flevel, its energy
*	as well as its EM-transition probability and its conversion coefficient.
*/
typedef struct gammatransition
{
	int ilevel;
	double energy;
	double probability;
	double cecoef;
	int flevel;
}GTrans;

double endpointenergy(double gammaenergy, double neutronenergy, int *w_sign);
void loadBetaFile();
void loadChargeFile();
void loadIsotopeFile(char* IsotopeFile);
void loadFermi(char* IsotopeFile);
double FF021p10(double p);
void loadEnergyFile(char *filename, ELevel elevels[], GTrans gammas[]);
int cascade(long int *seed, int startlevel,ELevel elevels[], GTrans gammas[], int *cascNum, double energies[],int *ceNum, double ce[]);