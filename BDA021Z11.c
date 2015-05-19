#include <stdio.h>
#include <math.h>
#include "BD2.h"
#include "BDA021Z11.h"
#include "c3ran.h"
#include "GaussianDev.h"

/*
 double *TOF		(output) ion time of flight for event
 double pe[4]	(output) positron 4-momentum in units of electron rest mass = 0.511 MeV
 double pv[4]	(output) neutrino 4-momentum in units of electron rest mass = 0.511 MeV
 double p1[4]	(output) 14O 4-momentum in units of electron rest mass = 0.511 MeV
 double pg[4]	(output) 2.3 MeV gamma-ray 4-momentum in units of electron rest mass
 double uran		(input) uniform random number generator function (I use Numerical Recipies ran2)
 long *sptr		(input) pointer to random number seed
 
 4-momentum:		p[0] = total energy = kinetic energy + rest mass
 p[1] = x component of momentum
 p[2] = y component of momentum
 p[3] = z component of momentum
 (rest mass) = sqrt(p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3])
 */
double trap_radius(double uran(long *),long *sptr);
double betadetector_hit(double *pe, double *loc, double *pe_exp, double uran(long *), long *sptr);
double MCPDetector_hit(double MRecoil, double *p1, double *loc, double *loc_a, double *TOF, double uran(long *), long *sptr);
double MCPDetector_accept(double E_n);
double VSCorr(double E,double E0);
double VSCorr2(double E,double E0);
double HCorr1(double E,double E0,double K,double pe_dot_k,int *w_sign);
double HCorr_a(double E,double E0,double K,double pe_dot_k,double pv_dot_k,double pe_dot_pv,int *w_sign);
double dGamma1(double E,double cos,double E0,double M,double cpex,double cpey,double cpez,double cpvx,
               double cpvy,double cpvz,int *w_sign);
double dGamma_a(double E,double cos,double E0,double M,double cpex,double cpey,double cpez,double cpvx,
                double cpvy,double cpvz,int *w_sign);
double RO1(double E,double E0,double M);
double ROa(double E,double E0,double M);
double EM1(double E,double E0,double M);
double EMa(double E,double E0,double M);
double dGammaH(double E,double E0,int *w_sign);

void BDA021Z11(double *TOF,int *w_sign,double *pe,double *pv,double *p1,double *pg,double pgarray[][4],double pcearray[][4], double *pn, 
               double *pe_exp,double *pv_exp,double *p1_exp,double *pg_exp,double *pn_exp,
               double *loc,double *loc_n,
               double *beta_rad, double *MCP_hit,double uran(long *),long *sptr,double M0exp,double M1exp,double M2exp,double Sn_exp,double abv,
			   int numGammas,double EGammas, double Gammas[],int cenum, double *ce, double Neutron1, int choice,double maxdist)
{
	static double m;			/* electron rest mass (MeV) */
	static double u;			/* conversion from amu to MeV */
	static double M0;			/* Iodine-137 rest mass (electron rest mass = 0.511 MeV) */
	static double M1;		        /* Xenon-137 rest mass (electron rest mass) */
	static double M2;		        /* Xenon-136 rest mass (electron rest mass) */
    static double m_n;                      /* neutron rest mass (in electron rest mass units) */
	static double Sn;                       /* neutron separation energy of Xenon-137 */
	static double E_n, KE_n;                /* neutron energy, kinetic energy */
	static double Egamma;                   /* gamma-ray energy */
	static double Egamma2;                   /* gamma-ray energy */
	static double pshake;                   /* shake-off electron energy */     
	static double E0main;			/* maximum beta total energy (electron rest mass = 0.511 MeV) */
	static double E0side;
	static double ft;			/* 14O ft-value (seconds) */
	static double w0main;			/* weight value (1/seconds) */
	static double w0side;
	static double maxEe_distmain;
	static double maxEe_distside;
	static int init=1;			/* initialization flag */
	
	double E0,w0,p,cex,cey,cez,cvx,cvy,cvz,cgx,cgy,cgz,csx,csy,csz,cnx,cny,cnz,ps[4];
	double rho0VS,rhoH;
	double costhetae,sinthetae,cosphie,sinphie;
	double cosbrem,sinbrem,phibrem,cosphibrem,sinphibrem,cbx,cby,cbz,pb[4],omega,beta,N;
	double pe_dot_k,pv_dot_k,pe_dot_pv;
	double dec_radius,a_betanu;
	double time,locb[3],cone_rad,rodchecky,rodcheckz,delta_rc;
	double dEe,Ee,Ee_dist,maxEe_dist,height,cosev,K1vsKa,cos_beta_neutron;
	double v_n,p_n,bn_height,beta_neutron_max;
    int i;
	static double a;
	FILE *rhofile;
    
	if (init==1) {
        // construct masses and energies in units of electron mass
		m=mexp;
		u=uexp;
		M0=1.0e-6*u*M0exp/m;
		M1=1.0e-6*u*M1exp/m;
		M2=1.0e-6*u*M2exp/m;
		m_n=1.0e-6*u*M_n_exp/m;
		//E_n=m_n+E_n_exp/m;
		Sn = Sn_exp/m;
		ft=ftexp;

		init = 0;
	}
		
		E_n=m_n+Neutron1/m;
		KE_n=E_n-m_n;
		a=abv;
		Egamma=EGammas/m;
		//printf("Gamma1: %f\n",Gamma1);
		Egamma2=Egammaexp2/m;
		pshake=sqrt(2.0*M0*Eshakeexp/m);
		E0main=M0-M1+1.0-Sn-KE_n-KE_n*m_n/M2-Egamma-Egamma2;
		//printf("E0  = %lf MeV\n",(E0main-1.0)*m);
		E0side=M0-M1-KE_n-KE_n*m_n/M2-Egamma-Egamma2;
		w0main=(E0main-1.0)*log(2.0)/ft;
		w0side=(E0side-1.0)*log(2.0)/ft;		

		/* this for loop finds the maximum of the Ee curve and integrates the phase space */
		Ee=1.0;
		rho0VS=0.0;
		maxEe_distmain=maxdist;
		/*dEe=(E0main-1.0)/(double)DIVISIONS;
        for(i=0;i<(DIVISIONS-1);i++) {
			Ee+=dEe;
            Ee_dist=dGamma1(Ee,1.0,E0main,M0,1.0,1.0,1.0,1.0,1.0,1.0,w_sign);
            if(Ee_dist > maxEe_distmain)
                maxEe_distmain=Ee_dist;
			rho0VS+=Ee_dist*dEe;
		}
		Ee=1.0;
		dEe=(E0side-1.0)/(double)DIVISIONS;
		for(i=0;i<(DIVISIONS-1);i++) {
			Ee+=dEe;
			Ee_dist=dGamma1(Ee,1.0,E0side,M0,1.0,1.0,1.0,1.0,1.0,1.0,w_sign);
			if(Ee_dist > maxEe_distside)
				maxEe_distside=Ee_dist;
        }
		Ee=1.0;
		rhoH=0.0;
		dEe=(E0main-1.0)/(double)DIVISIONS;
		for(i=0;i<(DIVISIONS-1);i++) {
			Ee+=dEe;
			Ee_dist=dGammaH(Ee,E0main,w_sign);
			rhoH+=Ee_dist*dEe;
		}
        
		//printf("En  = %lf MeV\nEr  = %lf MeV\nEg  = %lf MeV\nEg2 = %lf MeV\n",(KE_n)*m,(KE_n*m_n/M2)*m,Egamma*m,Egamma2*m);
		//printf("\n%lf	%lf\n\n",rho0VS,rhoH);
        //		rhofile=fopen("test.txt","w");
		*/
	if(choice==1 || choice==2) {
		E0=E0main;
		w0=w0main;
		a_betanu=a;
		maxEe_dist=maxEe_distmain;
		//if(choice==2)
			//maxEe_dist=BREMMAX;
	}
	if(choice==3 || choice==4) {
		E0=E0side;
		w0=w0side;
		a_betanu=(-1.0)*a_side;
		maxEe_dist=maxEe_distside;
		if(choice==4)
			maxEe_dist=BREMMAX;
	}
	if(choice<1 || choice>4) {
		printf("A branch error has occured\n\n");
		w0=-1.0;	/* this will remove the event */
		return;
	}
    
	for(;;) {
        /* generate uniform, random positron energy between 1 and E0 */
        pe[0]=(E0-1.0)*uran(sptr)+1.0;
        p=sqrt(pe[0]*pe[0]-1.0);
        
        /* generate uniform positron direction cosines */
        c3ran(&cex,&cey,&cez,uran,sptr);
        
        /* generate uniform neutrino direction cosines */
        c3ran(&cvx,&cvy,&cvz,uran,sptr);
		cosev = cex*cvx + cey*cvy + cez*cvz;
        
        /********************************************************************************************************/
		/* hard brem section -- see paper by F. Gluck in Computer Physics Communications 
         101 (1997) 223-231.  I tried to keep some of the same variables 
         for familiarity
         */
		/*
		if(choice==2 || choice==4) {
			omega=CS*(E0-pe[0]);
			beta=p/pe[0];
			N=0.5*log((1.0+beta)/(1.0-beta));
			pb[0]=omega*exp(-log(CS)*uran(sptr));
			cosbrem=(1.0-(1.0+beta)*exp(-2.0*N*uran(sptr)))/beta;
			sinbrem=sqrt(1.0-cosbrem*cosbrem);
			phibrem=2.0*PI*uran(sptr);
			cosphibrem=cos(phibrem);
			sinphibrem=sin(phibrem);
			pv[0]=E0-pe[0]-pb[0];
            
			costhetae=cez;
			sinthetae=sqrt(1.0-costhetae*costhetae);
			cosphie=cex/sinthetae;
			sinphie=cey/sinthetae;
            
			cbx=cex*cosbrem-sinphie*sinbrem*cosphibrem-costhetae*cosphie*sinphibrem*sinbrem;
			cby=cey*cosbrem+cosphie*sinbrem*cosphibrem-costhetae*sinphie*sinphibrem*sinbrem;
			cbz=cez*cosbrem+sinthetae*sinphibrem*sinbrem;
            
			pe_dot_k=p*pb[0]*cosbrem;
			pe_dot_pv=p*pv[0]*cosev;
			pv_dot_k=pv[0]*pb[0]*(cvx*cbx+cvy*cby+cvz*cbz);
            
			height=10.0*maxEe_dist*uran(sptr);
			rhoH=HCorr1(pe[0],E0,pb[0],pe_dot_k,w_sign)*2.0*PI*(E0-1.0)*log(1.0/CS);
			if(height < HCorr1(pe[0],E0,pb[0],pe_dot_k,w_sign)
			   + HCorr_a(pe[0],E0,pb[0],pe_dot_k,pv_dot_k,pe_dot_pv,w_sign)) {
                break;
			}
        */    
			/* end of hard brem section */
            /************************************************************************************************************/
		
		//} else {
            /* set hard bremsstrahlung photon energy to zero */
            pb[0]=0.0;
            cbx=0.0;
            cby=0.0;
            cbz=0.0;
            /* solve for neutrino energy by conservation of 4-momentum */
			pv[0]=(E0-pe[0]-(E0*E0-1.0)/(2.0*M0))/(1.0-(pe[0]-p*cosev)/M0);
			if(1) {
				/* generate random location between 0 and 2*maxEe_dist, where the 2 takes into account the 
                 contribution from the correlation term which can be just as big */
				height=2.0*maxEe_dist*uran(sptr);
				/* now we need to check whether this is under the graph -- if it is we 
                 will keep the point and this will give our positron energy -- if
                 not we will toss it out and generate another position */
				if(height < dGamma1(pe[0],cosev,E0,M0,cex,cey,cez,cvx,cvy,cvz,w_sign)
				   + a*dGamma_a(pe[0],cosev,E0,M0,cex,cey,cez,cvx,cvy,cvz,w_sign)) {
                    break;
				} 
			}
		//}
	}
    
	/* generate decay location in trap */
	loc[0]=trap_x+trap_radius(uran,sptr)*TRAP_FWHM_XSCALE;
	loc[1]=trap_y+trap_radius(uran,sptr)*TRAP_FWHM_YSCALE;
	loc[2]=trap_z+trap_radius(uran,sptr)*TRAP_FWHM_ZSCALE;
    
	/* electron momentum */
	pe[1]=p*cex;
	pe[2]=p*cey;
	pe[3]=p*cez;
    
	// does the beta hit the beta-detector?
	*beta_rad=betadetector_hit(pe,loc,pe_exp,uran,sptr);
    
	/* neutrino momentum */
	pv[1]=pv[0]*cvx;
	pv[2]=pv[0]*cvy;
	pv[3]=pv[0]*cvz;
	/* hard bremsstrahlung momentum */
    pb[1]=pb[0]*cbx;
    pb[2]=pb[0]*cby;
	pb[3]=pb[0]*cbz;
	/* shake-off electron momentum */
	c3ran(&csx,&csy,&csz,uran,sptr);
	ps[1]=pshake*csx;
	ps[2]=pshake*csy;
	ps[3]=pshake*csz;
	// if(choice==3 || choice==4) {
    /* generate random gamma direction */
	pg[0]=0;
	pg[1]=0;
	pg[2]=0;
	pg[3]=0;
	for(int j=0;j<numGammas;j++)
	{
		c3ran(&cgx,&cgy,&cgz,uran,sptr);
		pgarray[j][1]=Gammas[j]*cgx/m;
		pgarray[j][2]=Gammas[j]*cgy/m;
		pgarray[j][3]=Gammas[j]*cgz/m;
		pgarray[j][0]=sqrt(pgarray[j][1]*pgarray[j][1]+pgarray[j][2]*pgarray[j][2]+pgarray[j][3]*pgarray[j][3]);
		pg[1]+=pgarray[j][1];
		pg[2]+=pgarray[j][2];
		pg[3]+=pgarray[j][3];
		pg[0]+=pgarray[j][0];
	}
	// Conversion electron momentum generation
	// added 5/2 ksiegl
	for(int k=0;k<cenum;k++)
	{
		double pce=sqrt((ce[k]/m+1)*(ce[k]/m+1)-1);
		c3ran(&cgx,&cgy,&cgz,uran,sptr);
		pcearray[k][1]=pce*cgx;
		pcearray[k][2]=pce*cgy;
		pcearray[k][3]=pce*cgz;
		pcearray[k][0]=ce[k]/m;
		pg[1]+=pcearray[k][1];
		pg[2]+=pcearray[k][2];
		pg[3]+=pcearray[k][3];
		pg[0]+=pcearray[k][0];
	}
	/* borrowing shake-off for now */
	c3ran(&csx,&csy,&csz,uran,sptr);
	ps[1]=Egamma2*csx;
	ps[2]=Egamma2*csy;
	ps[3]=Egamma2*csz;
	ps[0]=sqrt(ps[1]*ps[1]+ps[2]*ps[2]+ps[3]*ps[3]);
	//} else {
	//  pg[0]=0.0;
	//  pg[1]=0.0;
	//  pg[2]=0.0;
	//  pg[3]=0.0;
	// }
	/* now calculate the recoil from all these contributions */
	p1[1]=(-(pe[1]+pv[1]+pb[1]+ps[1]+pg[1]));
	p1[2]=(-(pe[2]+pv[2]+pb[2]+ps[2]+pg[2]));
	p1[3]=(-(pe[3]+pv[3]+pb[3]+ps[3]+pg[3]));
	p1[0]=sqrt(p1[1]*p1[1]+p1[2]*p1[2]+p1[3]*p1[3] + M1*M1);
    
	/* now to calculate how recoil energy and momentum are affected */
    
	while(1) {
        /* generate random neutron direction */
        c3ran(&cnx,&cny,&cnz,uran,sptr);
        cos_beta_neutron=cex*cnx+cny*cny+cnz*cnz;
        beta_neutron_max=1.0+beta_neutron_p*pe[0]/(2.0*M1);        /* for possible beta-neutron angular correlation */
        bn_height=beta_neutron_max*uran(sptr);
        if(bn_height <= 1.0+(p/pe[0])*(p/pe[0])*beta_neutron_p*pe[0]/(2.0*M1)*cos_beta_neutron*cos_beta_neutron) {
            break;
        }
	}
    
    /* neutron momentum -- non-relativistic is good approximation */
	p_n = sqrt(E_n*E_n-m_n*m_n);
	v_n = p_n/m_n;
	pn[1] = m_n*(v_n*cnx + p1[1]/M1);
	pn[2] = m_n*(v_n*cny + p1[2]/M1);
	pn[3] = m_n*(v_n*cnz + p1[3]/M1);
	pn[0]=sqrt(pn[1]*pn[1]+pn[2]*pn[2]+pn[3]*pn[3]+m_n*m_n);
	
	if(KE_n>0.0) {
        /* recoil ion */
        p1[1] = -m_n*v_n*cnx + p1[1]*M2/M1;
        p1[2] = -m_n*v_n*cny + p1[2]*M2/M1;
        p1[3] = -m_n*v_n*cnz + p1[3]*M2/M1;
        p1[0] = sqrt(p1[1]*p1[1]+p1[2]*p1[2]+p1[3]*p1[3]+M2*M2);
	}

	/* does recoil ion a hit detector? returns hit location */
	*MCP_hit = MCPDetector_hit(M2,p1,loc,loc_n,TOF,uran,sptr);	
	if(*MCP_hit>0.0) {  // determine neutron energy from ion time of flight
        p1_exp[0] = M2*(loc_n[0]*loc_n[0]+loc_n[1]*loc_n[1]+loc_n[2]*loc_n[2])/(*TOF)/(*TOF)/SpeedofLight/SpeedofLight/2.0+M2;
        pn_exp[0] = (p1_exp[0]-M2)*M2/m_n + m_n;
	}
    
	*TOF=1.0/RF_FREQ*uran(sptr);
    if(*beta_rad<0.0 || *MCP_hit<0.0) {
		/* if we miss the beta detector or the HPGe detector, we want to signal the program to generate another
         event by making TOF < 0.0 */
		*TOF=-1.0;
        pe[0]=0.0;
        pe[1]=0.0;
        pe[2]=0.0;
        pe[3]=0.0;
        pv[0]=0.0;
        pv[1]=0.0;
        pv[2]=0.0;
        pv[3]=0.0;
        p1[0]=0.0;
        p1[1]=0.0;
        p1[2]=0.0;
        p1[3]=0.0;
	}
	return;
}

double trap_radius(double uran(long *), long *sptr)
{
	double radius,sigma;
	int i;
	/* for a Gaussian the relationship is FWHM=2.354*sigma */
	sigma=TRAP_FWHM/2.354;
	radius=GaussianDev(uran,sptr)*sigma;
	return(radius);
}

double betadetector_hit(double *pe, double *loc, double *pe_exp, double uran(long *), long *sptr)
{
    // right now I am assuming there is only one double-sided silicon strip detector located in the +x direction
	double x,y,z, Norm;	/* y,z position of where it would hit bd */
	double EnergySigma;
	
	return(1.0);  // this is in here so that all beta are considered

	// Beta in +z direction, with its active area in X-Y plane
	if(pe[3] > 0) {
	    x = (loc[0]-CENTER_X)+(pe[1]/pe[3])*(CENTER_Z-loc[2]+BD_Z);
	    y = (loc[1]-CENTER_Y)+(pe[2]/pe[3])*(CENTER_Z-loc[2]+BD_Z);
	    z = CENTER_Z-BD_Z;
	    // if(y<BD_SIDE/2.0 && y>-BD_SIDE/2.0 && z<BD_SIDE/2.0 && z>-BD_SIDE/2.0 && (pe[0]-1.0)>=BETA_THRESHOLD) {
	    if(x*x+y*y<BD_RADIUS*BD_RADIUS && (pe[0]-1.0)>=BETA_THRESHOLD) {
			x = (floor(x/BD_SIDE*NUM_STRIPS)+0.5)*BD_SIDE/NUM_STRIPS + (uran(sptr)-0.5)*BD_SIDE/NUM_STRIPS;
			y = (floor(y/BD_SIDE*NUM_STRIPS)+0.5)*BD_SIDE/NUM_STRIPS + (uran(sptr)-0.5)*BD_SIDE/NUM_STRIPS;
			
			// for a Gaussian the relationship is FWHM=2.354*sigma
			EnergySigma=sqrt(pe[0]-1.0)*BETA_ENERGY_FWHM/2.354;
			pe_exp[0] = pe[0] + EnergySigma*GaussianDev(uran,sptr);
			
			pe_exp[1] = x - CENTER_X;
			pe_exp[2] = y - CENTER_Y;
			pe_exp[3] = z - CENTER_Z;
			
			Norm = 1.0/sqrt(pe_exp[1]*pe_exp[1] + pe_exp[2]*pe_exp[2] + pe_exp[3]*pe_exp[3]);
			pe_exp[1] *= sqrt(pe_exp[0]*pe_exp[0] - 1.0)*Norm;
			pe_exp[2] *= sqrt(pe_exp[0]*pe_exp[0] - 1.0)*Norm;
			pe_exp[3] *= sqrt(pe_exp[0]*pe_exp[0] - 1.0)*Norm;
			return(1.0);
	    }
	}
	/*
    // Beta detector in -x direction, with its active area in Y-Z plane
	if(pe[1] < 0) {
	    x = CENTER_X-BD_X;
	    y = (loc[1]-CENTER_Y)+(pe[2]/pe[1])*(CENTER_X-loc[0]-BD_X);
	    z = (loc[2]-CENTER_Z)+(pe[3]/pe[1])*(CENTER_X-loc[0]-BD_X);
	    // if(y<BD_SIDE/2.0 && y>-BD_SIDE/2.0 && z<BD_SIDE/2.0 && z>-BD_SIDE/2.0 && (pe[0]-1.0)>=BETA_THRESHOLD) {
	    if(y*y+z*z<BD_RADIUS*BD_RADIUS && (pe[0]-1.0)>=BETA_THRESHOLD) {
			y = (floor(y/BD_SIDE*NUM_STRIPS)+0.5)*BD_SIDE/NUM_STRIPS + (uran(sptr)-0.5)*BD_SIDE/NUM_STRIPS;
			z = (floor(z/BD_SIDE*NUM_STRIPS)+0.5)*BD_SIDE/NUM_STRIPS + (uran(sptr)-0.5)*BD_SIDE/NUM_STRIPS;
			
			// for a Gaussian the relationship is FWHM=2.354*sigma
			EnergySigma=sqrt(pe[0]-1.0)*BETA_ENERGY_FWHM/2.354;
			pe_exp[0] = pe[0] + EnergySigma*GaussianDev(uran,sptr);
			
			pe_exp[1] = x - CENTER_X;
			pe_exp[2] = y - CENTER_Y;
			pe_exp[3] = z - CENTER_Z;
			
			Norm = 1.0/sqrt(pe_exp[1]*pe_exp[1] + pe_exp[2]*pe_exp[2] + pe_exp[3]*pe_exp[3]);
			pe_exp[1] *= sqrt(pe_exp[0]*pe_exp[0] - 1.0)*Norm;
			pe_exp[2] *= sqrt(pe_exp[0]*pe_exp[0] - 1.0)*Norm;
			pe_exp[3] *= sqrt(pe_exp[0]*pe_exp[0] - 1.0)*Norm;
			return(1.0);
	    }
	}
    */
    
    // Beta detector in +y direction, with its active area in X-Z plane
     /*
     if(pe[2] > 0) {
     x = (loc[0]-CENTER_X)-(pe[1]/pe[2])*(CENTER_Y-loc[1]-BD_Y);
     y = CENTER_Y-BD_Y; // this line doesn't matter
     z = (loc[2]-CENTER_Z)-(pe[3]/pe[2])*(CENTER_Y-loc[1]-BD_Y);
     // if(y<BD_SIDE/2.0 && y>-BD_SIDE/2.0 && z<BD_SIDE/2.0 && z>-BD_SIDE/2.0 && (pe[0]-1.0)>=BETA_THRESHOLD) {
     if(x*x+z*z<BD_RADIUS*BD_RADIUS && (pe[0]-1.0)>=BETA_THRESHOLD) {
     y = (floor(y/BD_SIDE*NUM_STRIPS)+0.5)*BD_SIDE/NUM_STRIPS + (uran(sptr)-0.5)*BD_SIDE/NUM_STRIPS;
     z = (floor(z/BD_SIDE*NUM_STRIPS)+0.5)*BD_SIDE/NUM_STRIPS + (uran(sptr)-0.5)*BD_SIDE/NUM_STRIPS;
     
     // for a Gaussian the relationship is FWHM=2.354*sigma 
    EnergySigma=sqrt(pe[0]-1.0)*BETA_ENERGY_FWHM/2.354;
    pe_exp[0] = pe[0] + EnergySigma*GaussianDev(uran,sptr);
    
    pe_exp[1] = x - CENTER_X;
    pe_exp[2] = y - CENTER_Y;
    pe_exp[3] = z - CENTER_Z;
    
    Norm = 1.0/sqrt(pe_exp[1]*pe_exp[1] + pe_exp[2]*pe_exp[2] + pe_exp[3]*pe_exp[3]);
    pe_exp[1] *= sqrt(pe_exp[0]*pe_exp[0] - 1.0)*Norm;
    pe_exp[2] *= sqrt(pe_exp[0]*pe_exp[0] - 1.0)*Norm;
    pe_exp[3] *= sqrt(pe_exp[0]*pe_exp[0] - 1.0)*Norm;
    return(1.0);
    }
    }
    */

return(-1.0);
}

double MCPDetector_hit(double m_2, double *p2, double *loc, double *loc_2, double *TOF, double uran(long *), long *sptr)
{
    // right now I am assuming there is only one MCP detector located in the +y direction and has same dimensions and location as DSSD
    double vel,dist;
    
    *TOF=1.0;
    return(1.0);  // all recoils are considered
    
    if(p2[2] > 0) {
        loc_2[0] = (loc[0]-CENTER_X)+(p2[1]/p2[2])*(BD_Y-loc[1]);
        loc_2[1] = BD_Y;
        loc_2[2] = (loc[2]-CENTER_Z)+(p2[3]/p2[2])*(BD_Y-loc[1]);
        
        // if(1) {
        if(loc_2[0]*loc_2[0]+loc_2[2]*loc_2[2] < 22.04*22.04) {
            // if(loc_2[0]<BD_SIDE/2.0 && loc_2[0]>-BD_SIDE/2.0 && loc_2[2]<BD_SIDE/2.0 && loc_2[2]>-BD_SIDE/2.0) {
            vel = sqrt(2.0*(p2[0]-m_2)/m_2)*SpeedofLight;
            dist = sqrt((loc_2[0]-loc[0])*(loc_2[0]-loc[0]) + (loc_2[1]-loc[1])*(loc_2[1]-loc[1]) + (loc_2[2]-loc[2])*(loc_2[2]-loc[2]));
            *TOF =  dist/vel;
            //if((dist/vel) > 0.0000000001)
            //printf("Time of flight: %lf   %lf   %lf\n", dist, vel, dist/vel*1.0e9);
            return(1.0);
        }
    }
    return(-1.0);
    
}

double MCPDetector_accept(double E_n)
{
	double acceptance;
    //	later we can add some functional form to the energy acceptance here
	acceptance=1.0;
	return(acceptance);
}

double VSCorr(double E,double E0)
{
	int i;
	double correction,spence,beta,x,N,omega;
	beta=sqrt(E*E-1.0)/E;
	x=2.0*beta/(1.0+beta);
	N=0.5*log((1.0+beta)/(1.0-beta));
	omega=CS*(E0-E);
	spence=0.0;
	for(i=1;i<VSACCURACY;i++) { 
		spence-=pow(x,i)/((double)i*(double)i);
	}
	correction=1.5*log(1837.4)+2.0*(N/beta-1.0)*log(2.0*omega)+2.0*N/beta*(1.0-N)+2.0/beta*spence-3.0/8.0;
	correction*=(ALPHA/PI);
	return(correction);
}

double VSCorr2(double E,double E0)
{
    double correction,beta,N;
    beta=sqrt(E*E-1.0)/E;
    N=0.5*log((1.0+beta)/(1.0-beta));
	correction=N*(beta*beta-1.0)/beta;
	correction*=(ALPHA/PI);
	return(correction);
}

double HCorr1(double E,double E0,double K,double pe_dot_k,int *w_sign)
{
	/* this function is the hard brem. function multiplied by the appropriate weight */
	double correction,H0;
	double beta,N,pe_4dot_k;
	beta=sqrt(E*E-1.0)/E;
	pe_4dot_k=E*K-pe_dot_k;
	N=0.5*log((1.0+beta)/(1.0-beta));
	H0=FF021p10(beta*E)*K*pe_4dot_k*(E0-E-K)
    *(-1.0*(E+K)*(1.0/K/K+1.0/pe_4dot_k/pe_4dot_k-2.0*E/K/pe_4dot_k)
      +K/pe_4dot_k);
	correction=(E0-E-K)*N*H0;
	correction*=(ALPHA/(2.0*PI*PI));
	if(correction>=0)
		*w_sign=1;
	else {
		*w_sign=-1;
		correction*=-1;
	}
	return(correction);
}

double HCorr_a(double E,double E0,double K,double pe_dot_k,double pv_dot_k,double pe_dot_pv,int *w_sign)
{
	double correction,H1;
    double beta,N,pe_4dot_k;
    beta=sqrt(E*E-1.0)/E;
	pe_4dot_k=E*K-pe_dot_k;
    N=0.5*log((1.0+beta)/(1.0-beta));
    H1=FF021p10(beta*E)*K*pe_4dot_k
    *(pe_dot_pv*(-1.0*(1.0/K/K+1.0/pe_4dot_k/pe_4dot_k-2.0*E/K/pe_4dot_k)+1.0/pe_4dot_k)
      +pv_dot_k*((E+K)/K/pe_4dot_k-1.0/pe_4dot_k/pe_4dot_k));
    correction=(E0-E-K)*N*H1;
    correction*=ALPHA/(2.0*PI*PI);
	if(correction>=0)
		*w_sign=1;
	else {
		*w_sign=-1;
		// correction*=-1;
	}
    return(correction);
}
double dGamma1(double E,double cosev,double E0,double M0,double cpex,double cpey,double cpez,double cpvx,double 
               cpvy, double cpvz,int *w_sign)
{
	/* this decay rate term covers isotropic term and A */
    double dgamma,pe;
    pe=sqrt(E*E-1.0);
    
	// simplify for now
    
	dgamma=FF021p10(pe)*pe*E*(E0-E)*(E0-E);  //uncomment this line to go back to an allowed decay
	
	//dgamma=FF021p10(pe)*pe*E*(E0-E)*(E0-E)*(((E*E)-1.0)+((E0-E)*(E0-E)));
    
	//	dgamma=FF021p10(pe)*pe*E*(E0-E)*(E0-E)*(1.0+VSCorr(E,E0)+VSCorr2(E,E0)+RO1(E,E0,M0)+EM1(E,E0,M0)
	//	+ bFierz/E + A*pe/E*(cpex*polarization[0]+cpey*polarization[1]+cpez*polarization[2]));
    
	if(dgamma>=0)
		*w_sign=1;
	else {
		*w_sign=-1;
		// dgamma*=-1;
	}
    return(dgamma);
}

double dGamma_a(double E,double cosev,double E0,double M0,double cpex,double cpey,double cpez,double cpvx,double cpvy,
                double cpvz,int *w_sign)
{
	/* this decay rate covers a, B, c, and D */
	double dgamma,pe;
	pe=sqrt(E*E-1.0);
    
	// simplify for now
    
	dgamma=FF021p10(pe)*pe*pe*(E0-E)*(E0-E)*cosev;  //uncomment this line to go back to allowed decay
	
	//dgamma=FF021p10(pe)*pe*pe*(E0-E)*(E0-E)*cosev*(((E*E)-1.0)+((E0-E)*(E0-E)));
    
    //	dgamma=FF021p10(pe)*pe*pe*(E0-E)*(E0-E)*cosev*(1.0+VSCorr(E,E0)+ROa(E,E0,M0)+EMa(E,E0,M0)
    //	+ B/(a*cosev)*E/pe*(cpvx*polarization[0]+cpvy*polarization[1]+cpvz*polarization[2]) 
    //	+ c_align/(a*cosev)*ALIGN*(cosev/3.0-(cpex*alignment[0]+cpey*alignment[1]+cpez*alignment[2])
    //		*(cpvx*alignment[0]+cpvy*alignment[1]+cpvz*alignment[2]))
    //	+ Dtrv/(a*cosev)*(polarization[0]*(cpey*cpvz-cpez*cpvy)-polarization[1]*(cpex*cpvz-cpez*cpvx)
    //		+polarization[2]*(cpex*cpvy-cpey*cpvx)));
	if(dgamma>=0)
		*w_sign=1;
	else {
		*w_sign=-1;
		// dgamma*=-1;
	}
	return(dgamma);
}

double dGammaH(double E,double E0,int *w_sign)
{
	double dgamma,pe,beta,N;
	pe=sqrt(E*E-1.0);
	beta=pe/E;
	N=0.5*log((1.0+beta)/(1.0-beta));
	dgamma=FF021p10(pe)*pe*E*(E0-E)*(E0-E)*(2.0*(N/beta-1.0)*(log(1.0/CS)+(E0-E)/3.0/E-1.5)
                                            +N*(E0-E)*(E0-E)/(12.0*beta*E*E));
	dgamma*=(ALPHA/PI);
	return(dgamma);
}

double RO1(double E,double E0,double M)
{
	double ROcorr1;
	ROcorr1 = 2.0*E0/(3.0*M)*MGT*(-MGT+dSCC+bWM) 
    + 2.0*E/(3.0*M)*(3.0*MF*MF+5.0*MGT*MGT-2.0*MGT*bWM)
    - 1.0/(3.0*M*E)*MGT*(2.0*MGT-dSCC-2.0*bWM);
	ROcorr1=ROcorr1/(MF*MF+MGT*MGT); /* beta spectrum is normalized such that MF*MF+MGT*MGT=1 */
	return(ROcorr1);
}

double ROa(double E,double E0,double M)
{
    double ROcorra;
	ROcorra = 2.0*E0/(3.0*M)*MGT*(MGT-dSCC-bWM)
    - 4.0*E/(3.0*M)*MGT*(3.0*MGT-bWM);
	ROcorra=ROcorra/(MF*MF-MGT*MGT/3.0);
	return(ROcorra);
}

double EM1(double E,double E0,double M)
{
    double EMCorr1;
	EMCorr1=6.0*ALPHA*Z*NUCRAD/35.0*(
                                     (MF*MF-MGT*MGT/3.0)*E0
                                     + (8.0*MF*MF+28.0*MGT*MGT/3.0)*E
                                     + 3.0/E*(MF*MF+MGT*MGT));
	EMCorr1=EMCorr1/(MF*MF+MGT*MGT);
	return(EMCorr1);
}

double EMa(double E,double E0,double M)
{
	double EMCorra;
	EMCorra=6.0*ALPHA*Z*NUCRAD/35.0*(
                                     (MF*MF+MGT*MGT)*E0
                                     + (8.0*MF*MF-4.0*MGT*MGT)*E);
	EMCorra=EMCorra/(MF*MF-MGT*MGT/3.0);
	return(EMCorra);
}
