#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "BD2.h"
#include "BDA021Z11.h"
#include "ran2.h"

#include "UVPSpline.h"
#include "CompUVPASpline.h"
#include "UVPSval.h"

struct betalevel
{
	double intensity;
	double gamma;
	double Ee_maxdist;
	int levelIndex; // Left this in here from Kevin, if he wants to use it.
};

struct neutronlevel
{
	double energy;
	double intensity;
	double Ee_maxdist;
};

struct chargestates
{
	int charge;
	double chargeweight;
};

struct isotope
{
	double M0;
	double M1;
	double M2;
	double Sn;
	double abv;
	double kbinding;
	int beta_num;
	int num_charge;
	int neutron_points;
	int foobar;
};

struct fermifunc
{
	double momentum;
	double fermi;
};


betalevel betaspec[beta_states];
neutronlevel neutronspec[neutron_states];
chargestates chargespec[charge_states];
isotope isoproperties;
fermifunc fermispec[fermi_points];
ELevel elevels[LVLNUM];
GTrans transitions[GTNUM];
double TRAP_FWHM=3.0;
int NCASE=100000; 	/* number of ions TOTAL! */
char* IsotopeFile=(char *)"134Sb_info.txt";
char* EnergyFile=(char *)"134TeRIPL.txt";

double pdat[fermi_points], Fdat[fermi_points];

int main(int argc, char *argv[])
{
	double TOF,pe[4],pv[4],p1[4],pg[4],pn[4];
	//added 4/1 ksiegl
	double pgarray[GAMMANUM][4];
	double pcearray[GAMMANUM][4];
	double Gammas[GAMMANUM];
	double ce[GAMMANUM];
	
	double pe_exp[4],pv_exp[4],p1_exp[4],pg_exp[4],pn_exp[4];
	double loc[3],M0,M1,M2,sumw1,sumwa,cosev,cosev_exp,beta_rad,MCP_hit;
	double loc_n[3];
	double me,m_n,pa1_p,pa2_p,Recoil,AlphaDetector[4],AlphaTOF[4];
	long i,j,seed;
	FILE *gamma_iondata, *gamma_geantdata, *gamma_auxdata, *neutron_iondata, *neutron_geantdata, *neutron_auxdata;
	int nscan;
	int choice,term,w_sign,num_files,charge,filetype,ntotal;
	char filenum[4],simion_filename[50],geant_filename[50];
	double Esum,pvsum,p1sum,p2sum,p3sum;
	double M1_SimIon, Theta_SimIon, Phi_SimIon, KE_SimIon;
	double M0exp, M1exp, M2exp, Sn_exp, abv;
	double Gamma1;
	double Ee_max;
	int state, charge_state, neutron_state;
	int gammanum;
	int cenum;
	double Egamma_weight_sum, neutron_weight_sum;
	double Eg_weight, neutron_weight;
	double chargeweight_sum, charge_weight;
	int decays, rundecays, eventnum, alldecays;
	float solidangle;
	int cetot,cestatetot;

	if(argc>1)
	{
		char temp[40];
		int length;
		IsotopeFile=argv[1];
		EnergyFile=argv[2];
		printf("load IsotopeFile %s\n",IsotopeFile);
		printf("load EnergyFile %s\n",EnergyFile);
		strcpy(temp,IsotopeFile);
		length = strlen(temp);
		if(temp[length-4]='.')
			temp[length-4]='\0';
		strcpy(GAMMA_SIMION_FILE,temp);
		strcat(GAMMA_SIMION_FILE,"_g_simion.ion");
		strcpy(GAMMA_GEANT_FILE,temp);
		strcat(GAMMA_GEANT_FILE,"_g_geant.txt");
		strcpy(GAMMA_AUXFILE,temp);
		strcat(GAMMA_AUXFILE,"_g_aux.txt");
		strcpy(NEUTRON_SIMION_FILE,temp);
		strcat(NEUTRON_SIMION_FILE,"_n_simion.ion");
		strcpy(NEUTRON_GEANT_FILE,temp);
		strcat(NEUTRON_GEANT_FILE,"_n_geant.txt");
		strcpy(NEUTRON_AUXFILE,temp);
		strcat(NEUTRON_AUXFILE,"_n_aux.txt");
	}
	
	if(argc>2)
	{
		NCASE=atoi(argv[3]);
	}
	
	if(argc>3)
	{
		choice=atoi(argv[4]);
	}
	else
	{
	printf("Press '1' for GAMMAS\n");
	printf("Press '2' for NEUTRONS\n");
	printf("Or hit '3' to exit this program.\n\n");

	nscan=scanf("%d",&choice);
	}
	
	if(argc>4)
	{
		TRAP_FWHM=atof(argv[5]);
	}
	
	switch(choice) {
    
	case 1:
		seed=-51988195;
		// Load the isotope masses and Sn into the struct

		loadIsotopeFile(IsotopeFile);
		
		// Load energy levels
		loadEnergyFile(EnergyFile,elevels,transitions);
		// Load Fermi tables

		loadFermi(IsotopeFile);
		
		cetot=0;
		cestatetot=0;
		for (int k=0; k<fermi_points; k++){
			pdat[k] = fermispec[k].momentum;
			Fdat[k] = fermispec[k].fermi;
		}

		M0=1.0e-6*uexp*isoproperties.M0/mexp;      /* M0 in units of mexp */
		M1=1.0e-6*uexp*isoproperties.M1/mexp;      /* M1 in units of mexp */
		M2=1.0e-6*uexp*isoproperties.M2/mexp;      /* M2 in units of mexp */

		//printf("m0: %f, m1:%f, m2:%f\n",isoproperties.M0,isoproperties.M1,isoproperties.M2);

		eventnum = 1;
		gammanum = (isoproperties.beta_num)-1; // This is a dummy variable now, but will be used in the future to keep track of how many gammas emitted per event

		//  Open the files to write appended SimIon and GEANT input data

		gamma_iondata = fopen(GAMMA_SIMION_FILE,"w");

		gamma_geantdata = fopen(GAMMA_GEANT_FILE,"w");

		gamma_auxdata = fopen(GAMMA_AUXFILE,"w");
		fprintf(gamma_auxdata,"eventnum, tob, num_state, energy_state, charge_state, decay_x, decay_y, decay_z\n"
								"nu_E, nu_px, nu_py, nu_pz\n"
								"beta_E, beta_px, beta_py, beta_pz\n"
								"rec_E, rec_px, rec_py, rec_pz\n"
								"gammanum\n"
								"gN_E, gN_px, gN_py, gN_pz\n"
								"CEnum\n"
								"ceN_E, ceN_px, ceN_py, ceN_pz\n"
								"neutron_num\n"
								"n_E, n_px, n_py, n_pz\n");

		// Get the sum of the beta intensities for normalization of decays to generate
		
		for(state=0;state<isoproperties.beta_num;state++){

			Egamma_weight_sum += betaspec[state].intensity;

		}

		//  Get the sum of the charge state weightings

		for(charge_state=0;charge_state<isoproperties.num_charge;charge_state++){

			chargeweight_sum += chargespec[charge_state].chargeweight;

		}

		//  Get the sum of the neutron intensities

		for(neutron_state=0;neutron_state<isoproperties.neutron_points;neutron_state++){

			neutron_weight_sum += neutronspec[neutron_state].intensity;
		
		}

		rundecays = 0;
		alldecays = 0;
		printf("beta states %i\n",isoproperties.beta_num);
		for(state=0;state<isoproperties.beta_num;state++){
			
			me=mexp;
			m_n=M_n_exp*1.0e-6*uexp/mexp;

			M0exp = isoproperties.M0;
			M1exp = isoproperties.M1;
			M2exp = isoproperties.M2;
			//Sn_exp = isoproperties.Sn;
			Sn_exp = 0.0;
			isoproperties.Sn = 0.0;
			abv = isoproperties.abv;

			Gamma1 = betaspec[state].gamma;

			decays = 0;

			//  Calculate the beta intensity weighting for each state, then generate the appropriate # of events

			Eg_weight = betaspec[state].intensity/Egamma_weight_sum;

			// Calculate the Q value for each state

			betaspec[state].Ee_maxdist=endpointenergy(Gamma1,0.0,&w_sign);

			for(charge_state=0;charge_state<isoproperties.num_charge;charge_state++){

				charge_weight = chargespec[charge_state].chargeweight/chargeweight_sum;
				
				ntotal = (int)ceil(NCASE*Eg_weight*charge_weight);

				printf("Eg_weight: %f, charge: %i, charge_weight: %f, ntotal: %i\n",Eg_weight,chargespec[charge_state].charge,charge_weight,ntotal);

				for(i=0;i<ntotal;++i) {
					// generate one decay -- not hitting both detectors means you have to generate yet another decay 
					for(j=0;j<4;j++) {
						pe[j] = pv[j] = p1[j] = pg[j] = pn[j] = 0.0;
						pe_exp[j] = pv_exp[j] = p1_exp[j] = pg_exp[j] = pn_exp[j] = 0.0;
					}

					Ee_max=betaspec[state].Ee_maxdist;
					
					//added 4/1 ksiegl
					gammanum=0.;
					cenum=0.;
					cascade(&seed, betaspec[state].levelIndex,elevels,transitions,&gammanum,Gammas,&cenum,ce);
					cestatetot+=cenum;
					//changed 4/1 ksiegl
					//printf("Cascade #%i energy %f\n",gammanum,betaspec[state].gamma);
					BDA021Z11(&TOF,&w_sign,pe,pv,p1,pg,pgarray,pcearray,pn,pe_exp,pv_exp,p1_exp,pg_exp,pn_exp,loc,loc_n,&beta_rad,&MCP_hit,ran2,&seed,
						M0exp,M1exp,M2exp,Sn_exp,abv,gammanum,Gamma1,Gammas,cenum,ce,0.0,choice,Ee_max);
					alldecays++;
					/* if particle doesn't hit detector TOF<0.0 */
					if(TOF>0.0) {
	
						Esum = M0-M1-(pe[0]-1.0)-(p1[0]-M1)-pg[0];
						p1sum = -(pe[1]+p1[1]+pg[1]);
						p2sum = -(pe[2]+p1[2]+pg[2]);
						p3sum = -(pe[3]+p1[3]+pg[3]);
						pvsum = sqrt(p1sum*p1sum+p2sum*p2sum+p3sum*p3sum);
	
						M1_SimIon=isoproperties.M1*1.0e-6-(chargespec[charge_state].charge+10*cenum)*mexp/uexp;
						KE_SimIon=(p1[0]-M1)*mexp*1e6;
					
						if(pn[0]-m_n>0.00001) {
							M1_SimIon=isoproperties.M2*1.0e-6-(chargespec[charge_state].charge+10*cenum)*mexp/uexp;
							KE_SimIon=(p1[0]-M2)*mexp*1e6;
						}
					
						Phi_SimIon=(-1.0)*atan(p1[3]/p1[1])*180.0/((double)PI);
						/* atan only returns values between 90 and -90 and SimIon only takes 
						values between 180 and -180.  Thus we must fix the ones that would 
						have larger than these values 
						*/
						if(p1[1]<0.0000000)
							Phi_SimIon+=180.0;
						if(Phi_SimIon>180.0)
							Phi_SimIon-=360.0;
						if(Phi_SimIon<-180.0)
							Phi_SimIon+=360.0;
	
						Theta_SimIon=atan(p1[2]/sqrt(p1[1]*p1[1]+p1[3]*p1[3]))*180.0/((double)PI);
					
						// Write to the SimIon and GEANT input files with all excites states appended
					
						fprintf(gamma_iondata, "%lf, ",TOF);
						fprintf(gamma_iondata, "%lf, ",M1_SimIon);
						fprintf(gamma_iondata, "%d, ",chargespec[charge_state].charge+10*cenum);
						fprintf(gamma_iondata, "%lf, ",loc[0]);
						fprintf(gamma_iondata, "%lf, ",loc[1]);
						fprintf(gamma_iondata, "%lf, ",loc[2]);
						fprintf(gamma_iondata, "%lf, ",Phi_SimIon);
						fprintf(gamma_iondata, "%lf, ",Theta_SimIon);
						fprintf(gamma_iondata, "%lf, ",KE_SimIon);
						fprintf(gamma_iondata, "1, ");
						fprintf(gamma_iondata, "1\n");
					
						fprintf(gamma_geantdata,"# %i\n",eventnum);
						fprintf(gamma_geantdata,"/gun/particle e-\n");
						fprintf(gamma_geantdata,"/gun/energy %lf MeV\n",(pe[0]-1.0)*mexp);
						fprintf(gamma_geantdata,"/gun/direction %lf %lf %lf\n",pe[1],pe[2],pe[3]);
						fprintf(gamma_geantdata,"/gun/position %lf %lf %lf mm\n",loc[0],loc[1],loc[2]);
						fprintf(gamma_geantdata,"/run/beamOn 1\n");

						// If there is a gamma ray emitted, write that information to the GEANT input file
						//changed 4/1 ksiegl
						if(gammanum>0){
							for(int k=0;k<gammanum;k++)
							{
								fprintf(gamma_geantdata,"/gun/particle gamma\n");
								fprintf(gamma_geantdata,"/gun/energy %lf MeV\n",Gammas[k]);
								fprintf(gamma_geantdata,"/gun/direction %lf %lf %lf\n",pgarray[k][1],pgarray[k][2],pgarray[k][3]);
								fprintf(gamma_geantdata,"/gun/position %lf %lf %lf mm\n",loc[0],loc[1],loc[2]);
								fprintf(gamma_geantdata,"/run/beamOn 1\n");
							}
						}

						// If there is a conversion electron emitted, write that information to the GEANT input file
						//added 5/2 ksiegl
						if(cenum>0){
							for(int k=0;k<cenum;k++)
							{
								fprintf(gamma_geantdata,"/gun/particle e-\n");
								fprintf(gamma_geantdata,"/gun/energy %lf MeV\n",ce[k]);
								fprintf(gamma_geantdata,"/gun/direction %lf %lf %lf\n",pcearray[k][1],pcearray[k][2],pcearray[k][3]);
								fprintf(gamma_geantdata,"/gun/position %lf %lf %lf mm\n",loc[0],loc[1],loc[2]);
								fprintf(gamma_geantdata,"/run/beamOn 1\n");
							}
						}
						
						// Write the auxiliary data file for post-analysis in ROOT
						// If running gamma rays and Sn = 0, set the neutron energy and momentum to 0 in output file

						if(Sn_exp == 0.0){
							pn[0] = 0.0;
							pn[1] = 0.0;
							pn[2] = 0.0;
							pn[3] = 0.0;
						}
						//changed 5/2 ksiegl
						fprintf(gamma_auxdata,"%i, %f, %i, %f, %i, %f, %f, %f\n%f, %f, %f, %f\n%f, %f, %f, %f\n%f, %f, %f, %f\n", 
								eventnum, TOF, state+1, betaspec[state].gamma, chargespec[charge_state].charge+10*cenum, loc[0], loc[1], loc[2], 
								pv[0], pv[1], pv[2], pv[3], 
								pe[0], pe[1], pe[2], pe[3], 
								p1[0], p1[1], p1[2], p1[3]); 
						fprintf(gamma_auxdata,"%i\n",gammanum);
						for(int k=0;k<gammanum;k++)
						{
							fprintf(gamma_auxdata,"%f, %f, %f, %f\n",pgarray[k][0],pgarray[k][1],pgarray[k][2],pgarray[k][3]);
						}
						fprintf(gamma_auxdata,"%i\n",cenum);
						for(int k=0;k<cenum;k++)
						{
							fprintf(gamma_auxdata,"%f, %f, %f, %f\n",pcearray[k][0],pcearray[k][1],pcearray[k][2],pcearray[k][3]);
						}
						fprintf(gamma_auxdata,"0\n");
					
						eventnum++;
					}
					else if (TOF<0.0)
						i--;

				}

				decays += ntotal;
				cetot+=cestatetot;
				printf("%d conversion electrons\n",cestatetot);
				cestatetot=0;
			}
		
		rundecays += decays;
		printf("%d decays were generated for E_Ex = %f MeV.\n",decays,Gamma1);
		}
		fclose(gamma_iondata);
		fclose(gamma_geantdata);
		fclose(gamma_auxdata);
		printf("-----------------SUMMARY-----------------\n");
		printf("%d decays generated for this isotope.\n",rundecays);
		printf("with %d conversion electrons.\n",cetot);
		solidangle = (float(rundecays)/float(alldecays))*100;
		printf("%d decays generated to hit beta detector.\n",alldecays);
		printf("Beta detector has a %0.2f%% solid angle.\n",solidangle);
		break;
	case 2:
		seed=-51988195;

		// Load the isotope masses and Sn into the struct
		isoproperties.neutron_points=0;
		loadIsotopeFile(IsotopeFile);

		// Load Fermi tables

		loadFermi(IsotopeFile);

		for (int k=0; k<fermi_points; k++){
			pdat[k] = fermispec[k].momentum;
			Fdat[k] = fermispec[k].fermi;
		}

		M0=1.0e-6*uexp*isoproperties.M0/mexp;      /* M0 in units of mexp */
		M1=1.0e-6*uexp*isoproperties.M1/mexp;      /* M1 in units of mexp */
		M2=1.0e-6*uexp*isoproperties.M2/mexp;      /* M2 in units of mexp */

		//printf("m0: %f, m1:%f, m2:%f\n",isoproperties.M0,isoproperties.M1,isoproperties.M2);

		eventnum = 1;
		gammanum = 0; // This is a dummy variable now, but will be used in the future to keep track of how many gammas emitted per event

		//  Open the files to write appended SimIon and GEANT input data

		neutron_iondata = fopen(NEUTRON_SIMION_FILE,"w");

		neutron_geantdata = fopen(NEUTRON_GEANT_FILE,"w");

		neutron_auxdata = fopen(NEUTRON_AUXFILE,"w");
		fprintf(neutron_auxdata,"eventnum, tob, num_state, energy_state, charge_state, decay_x, decay_y, decay_z\n"
								"nu_E, nu_px, nu_py, nu_pz\n"
								"beta_E, beta_px, beta_py, beta_pz\n"
								"rec_E, rec_px, rec_py, rec_pz\n"
								"gammanum\n"
								"gN_E, gN_px, gN_py, gN_pz\n"
								"CEnum\n"
								"ceN_E, ceN_px, ceN_py, ceN_pz\n"
								"neutron_num\n"
								"n_E, n_px, n_py, n_pz\n");
		//fprintf(neutron_auxdata,"eventnum, tob, num_state, energy_state, charge_state, decay_x, decay_y, decay_z, nu_E, nu_px, nu_py, nu_pz, beta_E, beta_px, beta_py, beta_pz, rec_E, rec_px, rec_py, rec_pz, gammanum, g1_E, g1_px, g1_py, g1_pz, n_E, n_px, n_py, n_pz\n");

		//  Get the sum of the neutron intensities

		for(neutron_state=0;neutron_state<isoproperties.neutron_points;neutron_state++){

			neutron_weight_sum += neutronspec[neutron_state].intensity;
		
		}

		//  Get the sum of the charge state weightings
		
		for(charge_state=0;charge_state<isoproperties.num_charge;charge_state++){

			chargeweight_sum += chargespec[charge_state].chargeweight;

		}
		
		rundecays = 0;
		alldecays = 0;
		printf("neutron states: %d\n",isoproperties.neutron_points);
		for(neutron_state=0;neutron_state<isoproperties.neutron_points;neutron_state++){
			
			me=mexp;
			m_n=M_n_exp*1.0e-6*uexp/mexp;

			M0exp = isoproperties.M0;
			M1exp = isoproperties.M1;
			M2exp = isoproperties.M2;
			Sn_exp = isoproperties.Sn;
			abv = isoproperties.abv;

			Gamma1 = 0.0; // For generating decays to neutron states, gamma energy needs to be 0.
			gammanum = 0;
			decays = 0;

			//  Calculate the neutron intensity weighting to generate the appropriate # of events

			neutron_weight = neutronspec[neutron_state].intensity/neutron_weight_sum;

			// Calculate the Q value for each state

			neutronspec[neutron_state].Ee_maxdist=endpointenergy(0,neutronspec[neutron_state].energy,&w_sign);
			
			for(charge_state=0;charge_state<isoproperties.num_charge;charge_state++){
				printf("Charge state: %i\n",chargespec[charge_state].charge);
				charge_weight = chargespec[charge_state].chargeweight/chargeweight_sum;
				
				ntotal = (int)ceil(NCASE*neutron_weight*charge_weight);

				printf("neutron_weight: %f, charge: %i, charge_weight: %f, ntotal: %i\n",neutron_weight,chargespec[charge_state].charge,charge_weight,ntotal);
				
				for(i=0;i<ntotal;++i) {
					// generate one decay -- not hitting both detectors means you have to generate yet another decay 
					for(j=0;j<4;j++) {
						pe[j] = pv[j] = p1[j] = pg[j] = pn[j] = 0.0;
						pe_exp[j] = pv_exp[j] = p1_exp[j] = pg_exp[j] = pn_exp[j] = 0.0;
					}

					Ee_max=neutronspec[neutron_state].Ee_maxdist;
				
					BDA021Z11(&TOF,&w_sign,pe,pv,p1,pg,pgarray,pcearray,pn,pe_exp,pv_exp,p1_exp,pg_exp,pn_exp,loc,loc_n,&beta_rad,&MCP_hit,ran2,&seed,
						M0exp,M1exp,M2exp,Sn_exp,abv,gammanum,Gamma1,Gammas,cenum,ce,neutronspec[neutron_state].energy,choice,Ee_max);
					
					/* if particle doesn't hit detector TOF<0.0 */
					if(TOF>0.0) {
	
						Esum = M0-M1-(pe[0]-1.0)-(p1[0]-M1)-pg[0];
						p1sum = -(pe[1]+p1[1]+pg[1]);
						p2sum = -(pe[2]+p1[2]+pg[2]);
						p3sum = -(pe[3]+p1[3]+pg[3]);
						pvsum = sqrt(p1sum*p1sum+p2sum*p2sum+p3sum*p3sum);
	
						M1_SimIon=isoproperties.M1*1.0e-6-chargespec[charge_state].charge*mexp/uexp;
						KE_SimIon=(p1[0]-M1)*mexp*1e6;
					
						if(pn[0]-m_n>0.00001) {
							M1_SimIon=isoproperties.M2*1.0e-6-chargespec[charge_state].charge*mexp/uexp;
							KE_SimIon=(p1[0]-M2)*mexp*1e6;
						}
					
						Phi_SimIon=(-1.0)*atan(p1[3]/p1[1])*180.0/((double)PI);
						/* atan only returns values between 90 and -90 and SimIon only takes 
						values between 180 and -180.  Thus we must fix the ones that would 
						have larger than these values 
						*/
						if(p1[1]<0.0000000)
							Phi_SimIon+=180.0;
						if(Phi_SimIon>180.0)
							Phi_SimIon-=360.0;
						if(Phi_SimIon<-180.0)
							Phi_SimIon+=360.0;
	
						Theta_SimIon=atan(p1[2]/sqrt(p1[1]*p1[1]+p1[3]*p1[3]))*180.0/((double)PI);
					
						// Write to the SimIon and GEANT input files with all excites states appended
					
						fprintf(neutron_iondata, "%lf, ",TOF);
						fprintf(neutron_iondata, "%lf, ",M1_SimIon);
						fprintf(neutron_iondata, "%d, ",chargespec[charge_state].charge);
						fprintf(neutron_iondata, "%lf, ",loc[0]);
						fprintf(neutron_iondata, "%lf, ",loc[1]);
						fprintf(neutron_iondata, "%lf, ",loc[2]);
						fprintf(neutron_iondata, "%lf, ",Phi_SimIon);
						fprintf(neutron_iondata, "%lf, ",Theta_SimIon);
						fprintf(neutron_iondata, "%lf, ",KE_SimIon);
						fprintf(neutron_iondata, "1, ");
						fprintf(neutron_iondata, "1\n");
					
						fprintf(neutron_geantdata,"# %i\n",eventnum);
						fprintf(neutron_geantdata,"/gun/particle e-\n");
						fprintf(neutron_geantdata,"/gun/energy %lf MeV\n",(pe[0]-1.0)*mexp);
						fprintf(neutron_geantdata,"/gun/direction %lf %lf %lf\n",pe[1],pe[2],pe[3]);
						fprintf(neutron_geantdata,"/gun/position %lf %lf %lf mm\n",loc[0],loc[1],loc[2]);
						fprintf(neutron_geantdata,"/run/beamOn 1\n");

						// If there is a gamma ray emitted, write that information to the GEANT input file

						if(neutronspec[neutron_state].energy>0){
							fprintf(neutron_geantdata,"/gun/particle neutron\n");
							fprintf(neutron_geantdata,"/gun/energy %lf MeV\n",neutronspec[neutron_state].energy);
							fprintf(neutron_geantdata,"/gun/direction %lf %lf %lf\n",pn[1],pn[2],pn[3]);
							fprintf(neutron_geantdata,"/gun/position %lf %lf %lf mm\n",loc[0],loc[1],loc[2]);
							fprintf(neutron_geantdata,"/run/beamOn 1\n");
						}

						// Write the auxiliary data file for post-analysis in ROOT

						//fprintf(neutron_auxdata,"%i, %f, %i, %f, %i, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %i,  %f, %f, %f, %f, %f, %f, %f, %f\n", eventnum, TOF, state+1, betaspec[state].gamma, chargespec[charge_state].charge, loc[0], loc[1], loc[2], pv[0], pv[1], pv[2], pv[3], pe[0], pe[1], pe[2], pe[3], p1[0], p1[1], p1[2], p1[3], gammanum, pg[0], pg[1], pg[2], pg[3], pn[0], pn[1], pn[2], pn[3]);
						
						fprintf(neutron_auxdata,"%i, %f, %i, %f, %i, %f, %f, %f\n%f, %f, %f, %f\n%f, %f, %f, %f\n%f, %f, %f, %f\n", 
								eventnum, TOF, state+1, betaspec[state].gamma, chargespec[charge_state].charge, loc[0], loc[1], loc[2], 
								pv[0], pv[1], pv[2], pv[3], 
								pe[0], pe[1], pe[2], pe[3], 
								p1[0], p1[1], p1[2], p1[3]); 
						fprintf(neutron_auxdata,"%i\n",0);
						fprintf(neutron_auxdata,"%i\n",0);
						fprintf(neutron_auxdata,"1\n");
						fprintf(neutron_auxdata,"%f, %f, %f, %f\n",pn[0],pn[1],pn[2],pn[3]);
						eventnum++;
					}
					else if (TOF<0.0)
						i--;

				}

				decays += ntotal;
			}
		
		rundecays += decays;
		printf("%d decays were generated for En = %f MeV.\n",decays,neutronspec[neutron_state].energy);
		}
		fclose(neutron_iondata);
		fclose(neutron_geantdata);
		fclose(neutron_auxdata);
		printf("-----------------SUMMARY-----------------\n");
		printf("%d decays generated for this isotope.\n",rundecays);
		solidangle = (float(rundecays)/float(alldecays))*100;
		printf("%d decays generated to hit beta detector.\n",alldecays);
		printf("Beta detector has a %0.2f%% solid angle.\n",solidangle);
		break;
	case 3:
		printf("Alright.  Take it easy then, dude.\n\n");
		break;
		default:
		break;
	}
	return(0);
}

double endpointenergy(double gammaenergy, double neutronenergy, int *w_sign)
{
	double m=mexp;
	double u=uexp;
	double M0=1.0e-6*u*isoproperties.M0/m;
	double M1=1.0e-6*u*isoproperties.M1/m;
	double M2=1.0e-6*u*isoproperties.M2/m;
	double m_n=1.0e-6*u*M_n_exp/m;
	//double E_n=m_n+E_n_exp/m;
	double E_n=m_n+neutronenergy/m;
	double KE_n=E_n-m_n;
	double Sn = isoproperties.Sn/m;
	double E0main=M0-M1+1.0-Sn-KE_n-KE_n*m_n/M2-gammaenergy/m;
	double E0side=M0-M1-KE_n-KE_n*m_n/M2-gammaenergy/m;
	double w0main=(E0main-1.0)*log(2.0)/ftexp;
	double w0side=(E0side-1.0)*log(2.0)/ftexp;		
	double Ee_dist=0;
	double maxEe_distmain=0;
		/* this for loop finds the maximum of the Ee curve and integrates the phase space */
	double Ee=1.0;
	double rho0VS=0.0;
	double dEe=(E0main-1.0)/(double)DIVISIONS;
    for(int i=0;i<(DIVISIONS-1);i++) {
		Ee+=dEe;
        Ee_dist=dGamma1(Ee,1.0,E0main,M0,1.0,1.0,1.0,1.0,1.0,1.0,w_sign);
        if(Ee_dist > maxEe_distmain)
			maxEe_distmain=Ee_dist;
		rho0VS+=Ee_dist*dEe;
	}
	//printf("m_n: %f\n",m_n);
	printf("E0  = %lf MeV\nEn  = %lf MeV\nEr  = %lf MeV\nEg  = %lf MeV\n",(E0main-1.0)*m,(KE_n)*m,(KE_n*m_n/M2)*m,gammaenergy);
	return maxEe_distmain;
}

void loadIsotopeFile(char* IsotopeFile)
{
	FILE * isotopefile;
	
	char s1[50], s2[50], s3[50], s4[50];
	char line[100];

	int num_charge, num_beta, num_neutron;

	double kbinding;

	isotopefile = fopen(IsotopeFile,"r");
	if (isotopefile == NULL) perror ("Error opening isotope.txt file");
	
	while(fgets(line,sizeof line, isotopefile) != NULL){
		//if (*line =='#') continue;
		if (sscanf(line, "M0: %lf, M1: %lf, M2: %lf, Sn: %lf", &isoproperties.M0, &isoproperties.M1, &isoproperties.M2, &isoproperties.Sn) == 4){
			printf("M0: %lf, M1: %lf, M2: %lf, Sn: %lf\n", isoproperties.M0, isoproperties.M1, isoproperties.M2, isoproperties.Sn);
		}
		if (sscanf(line, "abv: %lf", &isoproperties.abv) == 1){
			printf("abv: %lf\n", isoproperties.abv);
		}
		if (sscanf(line, "k-shell BE: %lf", &isoproperties.kbinding) == 1){
			printf("k-shell BE: %lf keV\n", isoproperties.kbinding);
		}
		if (sscanf(line, "Number of charge states: %i", &isoproperties.num_charge) == 1){
			printf("Number of charge states: %i\n", isoproperties.num_charge);
		}
		if (sscanf(line, "CHARGE\tWEIGHT%[^\n]", s1) == 1){
			//printf("CHARGE\tWEIGHT\n");
			fgets(line,sizeof line,isotopefile); // This fgets skips the CHARGE WEIGHT identifier
			for (int i=0;i<isoproperties.num_charge;i++){
				sscanf(line, "%i\t%lf", &chargespec[i].charge, &chargespec[i].chargeweight);
				printf("CHARGE: %i, WEIGHT: %lf\n", chargespec[i].charge, chargespec[i].chargeweight);
				fgets(line,sizeof line,isotopefile);
			}
		}

		if (sscanf(line, "Number of beta branches: %i", &num_beta) == 1){
			isoproperties.beta_num=num_beta;
			printf("Number of beta branches: %i\n", num_beta);
		}
		
		if (sscanf(line, "EXCITATION\tINTENSITY%[^\n]", s2) == 1){
			//printf("EXCITATION\tINTENSITY\n");
			fgets(line,sizeof line,isotopefile);
			for (int i=0;i<num_beta;i++){
				//changed 4/1 ksiegl
				sscanf(line, "%lf\t%lf\t%i", &betaspec[i].gamma, &betaspec[i].intensity,&betaspec[i].levelIndex);
				printf("EXCITATION: %lf, INTENSITY: %lf\n", betaspec[i].gamma, betaspec[i].intensity);
				fgets(line,sizeof line,isotopefile);
			}
		}

		if (sscanf(line, "Number of neutron states: %i", &num_neutron) == 1){
			printf("Number of neutron states: %i\n", num_neutron);
			isoproperties.neutron_points=num_neutron;
		}
		
		if (sscanf(line, "ENERGY\tINTENSITY%[^\n]", s3) == 1){
			fgets(line,sizeof line,isotopefile);
			for (int i=0;i<num_neutron;i++){
				sscanf(line, "%lf\t%lf", &neutronspec[i].energy, &neutronspec[i].intensity);
				printf("ENERGY: %lf, INTENSITY: %lf\n", neutronspec[i].energy, neutronspec[i].intensity);
				fgets(line,sizeof line,isotopefile);
			}
		}
		
	}
	isoproperties.neutron_points=num_neutron;
	printf("charge 0 %i\n",chargespec[0].charge);
	fclose(isotopefile);
	return;
}

void loadFermi(char* IsotopeFile)
{
	FILE * isotopefile;
	
	char s1[50];
	char line[100];

	isotopefile = fopen(IsotopeFile,"r");
	if (isotopefile == NULL) perror ("Error opening isotope.txt file");
	
	while(fgets(line,sizeof line, isotopefile) != NULL){

		if (sscanf(line, "MOMENTUM\tFermi%[^\n]", s1) == 1){
			//printf("MOMENTUM\tFermi\n");
			fgets(line,sizeof line,isotopefile);
			for (int i=0;i<fermi_points;i++){
				sscanf(line, "%lf\t%lf", &fermispec[i].momentum, &fermispec[i].fermi);
				//printf("MOMENTUM: %lf, Fermi: %lf\n", fermispec[i].momentum, fermispec[i].fermi);
				fgets(line,sizeof line,isotopefile);
			}
		}
		
	}
	
	fclose(isotopefile);
	return;
}

/*	Loads the RIPL-3 formatted energy level tables into a pair of arrays, one for
*	energy levels, the other for the associated gamma transitions.
*	
*	Currently it works, though there are some issues that need to be addressed regarding 
*	the presence of halflives for only some of the states.
*/
void loadEnergyFile(char *filename, energylevel elevels[], gammatransition gammas[])
{
	FILE *ef;
	char ignore[1024];
	char gammae[20];
	char gammap[20];
	char cecoef[20];
	int numLevel=0;
	int numGamma;
	int massNum;
	float levelEnergy;
	int levelGnum;
	float neutronSep;
	char check[10];
	int gammaIndex;
	int i;
	int flevel;
	double energy;
	double p;
	gammaIndex=0;
	printf("loading %s\n", filename);
	ef = fopen(filename,"r");
	printf("loaded\n");
	fscanf(ef,"%s",ignore);
	printf("Daughter Isotope %s\n", ignore);
	fscanf(ef,"%s %s %s %i",ignore,ignore, ignore, &numLevel);
	printf("Number %s %i\n", ignore, numLevel);
	fscanf(ef," number of gamma-rays: %i",&numGamma);
	printf("Number Transitions: %i\n",numGamma);
	fscanf(ef," neutron separation energy: %f %s",&neutronSep, ignore);
	//printf("fscanf check %s\n",ignore);
	for(i=0;i<9;i++) fgets(ignore,sizeof(ignore),ef);
	//printf("fscanf check ignore %s\n",ignore);
	for(int j=1;j<=numLevel;j++)
	{
		printf("level %i\n",j);
		int levelIndex;
		//Each level lists its index and its energy
		fscanf(ef,"%i %f %*f %*i",&levelIndex, &levelEnergy);
		fscanf(ef,"%s",check);
		//Because some levels have halflives listed in a "disappearing" column before transition number, a check for whether
		//the column is an integer or out of place (larger than the level number) is used.  This solution is imperfect.
		if(check[1]=='.'){
			fscanf(ef,"%i",&levelGnum);
		}
		else{
			levelGnum=atoi(check);
		}
		printf("Level gnum %i\n",levelGnum);
		elevels[levelIndex].numofgamma=levelGnum;
		elevels[levelIndex].energy=levelEnergy;
		elevels[levelIndex].firstgamma=gammaIndex;
		fgets(ignore,sizeof(ignore),ef);
		for(i=0;i<levelGnum;i++){
			fscanf(ef,"%i %s %*s %s %s",&flevel,gammae,gammap,cecoef);
			printf("fscanf check probability %s\n",gammap);
			gammas[gammaIndex].ilevel=levelIndex;
			gammas[gammaIndex].flevel=flevel;
			gammas[gammaIndex].energy=atof(gammae);
			gammas[gammaIndex].probability=atof(gammap);
			printf("Transition prob: %f\n",gammas[gammaIndex].probability);
			gammas[gammaIndex].cecoef=atof(cecoef);
			elevels[levelIndex].totprob=elevels[levelIndex].totprob+atof(gammap);
			gammaIndex++;
			fgets(ignore,sizeof(ignore),ef);
		}
	}	
}

/*	Randomly generates a gamma cascade from starting level startlevel in elevels[] with transitions gammas[].
*	The length of the generated cascade is passed through cascNum, the sequential transition energies are
*	stored in energies[] and whether it is an internal conversion is stored in ce[]
*/ 
int cascade(long int *seed, int startlevel,ELevel elevels[], GTrans gammas[], int *cascNum, double energies[],int *ceNum, double ce[])
{
	int level=startlevel;
	*cascNum=0;
	*ceNum=0;
	//printf("casc check %i %f\n",level,elevels[level].totprob);
	while(elevels[level].totprob>=1e-8)
	{
		//printf("level %i\n",level);
		float probnorm=elevels[level].totprob;
		float u=probnorm*ran2(seed);
		//printf("%i ran %f\n",level, u);
		int gamma1=elevels[level].firstgamma;
		int i=0;
		while(u>0. && i<elevels[level].numofgamma)
		{
			u=u-gammas[gamma1+i].probability;
			i++;
		}
		float icran=ran2(seed);
		float icc=gammas[gamma1+i-1].cecoef;
		//printf("energy %f\n",energies[*cascNum]);
		if(icran<(icc/(icc+1)))
		{
			ce[*ceNum]=gammas[gamma1+i-1].energy-isoproperties.kbinding/1000.;
			*ceNum = *ceNum+1;
			printf("conversion electron energy: %f\n",ce[(*ceNum)-1]);
		}
		else
		{
			energies[*cascNum]=gammas[gamma1+i-1].energy;
			*cascNum=*cascNum+1;
		}
		level=gammas[gamma1+i-1].flevel;
	}
	//printf("Cascade from %i to %i with %i transitions\n", startlevel, level, *cascNum);
	return 0;
}

static UVPSpline *sptr=NULL;

double FF021p10(double p)
{
	double Fval;
	/*
	double pdat[fermi_points], Fdat[fermi_points];

	for (int k=0; k<fermi_points; k++){
		pdat[k] = fermispec[k].momentum;
		Fdat[k] = fermispec[k].fermi;
	}
	*/
	if (!sptr) {
		sptr=sptr=CompUVPASpline(fermi_points,pdat,Fdat);
	}
	
	Fval=UVPSval(p,sptr);

	//printf("pdat: %lf, Fdat: %lf, Fval: %lf\n",p, Fdat, Fval);
	
	return (Fval>0.0)?Fval:0.0;
}