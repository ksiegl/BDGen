#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "BD2.h"
#include "GaussianDev.h"

void histogram(FILE *file,int counts[NUM_BE_BINS][NUM_GE_BINS][NUM_GE],
	double average[NUM_BE_BINS][2][NUM_GE],int filetype);
void histwrite(int K1[NUM_BE_BINS][NUM_GE_BINS][NUM_GE],int Ka[NUM_BE_BINS][NUM_GE_BINS][NUM_GE],
	double A1[NUM_BE_BINS][2][NUM_GE],double Aa[NUM_BE_BINS][2][NUM_GE]);
double Ge_res(double uran(long *),long *sptr);
double Si_res(double E,double uran(long *),long *sptr);
double gaussian(double E,double Eprime);
double lowEtail(double E,double Eprime);
double highEplateau(double E,double Eprime);
double highEtail(double E,double Eprime);

void BDSimFile(double uran(long *),long *sptr)
{
	FILE *input,*output,*output2;
	char return1,return2;
	char num[4],inputname[20],outputname[20],outputname2[20];
	int K1count[NUM_BE_BINS][NUM_GE_BINS][NUM_GE],Kacount[NUM_BE_BINS][NUM_GE_BINS][NUM_GE];
	double average_1[NUM_BE_BINS][2][NUM_GE],average_a[NUM_BE_BINS][2][NUM_GE];
	int i,j,k,nscan,nscan1,nscan2;
	int filenum,filetype,sign,choice,file1ora,ntotal;
	double pe0,pe1,pe2,pe3,pg0,pg1,pg2,pg3;
	double signvalue_1,sign_1,signvalue_a,sign_a;


	/* this file is not optimized for speed because it is so fast anyway --
	   to optimize the roles of i and j should be switched because of the 
	   way computers read in 2D arrays.  This would speed it up x2.
	*/

	for(i=0;i<NUM_BE_BINS;i++) {
		for(j=0;j<NUM_GE_BINS;j++) {
			for(k=0;k<NUM_GE;k++) {
				K1count[i][j][k]=0;
				Kacount[i][j][k]=0;
			}
		}
	}

	for(i=0;i<NUM_BE_BINS;i++) {
		for(j=0;j<2;j++) {
			for(k=0;k<NUM_GE;k++) {
				average_1[i][j][k]=0.0;
				average_a[i][j][k]=0.0;
			}
		}
	}

	for(filetype=0;filetype<NUM_GE;filetype++) {
		choice=1;
		ntotal=NCASE;
		if(choice==2)
		        ntotal=(int)floor(ntotal*HARDBR);
                if(choice==3)
			ntotal=(int)floor(ntotal*SIDEBR/(1.0-(double)SIDEBR));
		if(choice==4)
			ntotal=(int)floor(ntotal*SIDEBR/(1.0-(double)SIDEBR)*SIDEHARDBR);

		for(file1ora=0;file1ora<=1;file1ora++) {
			sprintf(num,"%d\0",filetype+1);
			strcpy(inputname,num);
			strcpy(outputname,num);

			if(file1ora==0) {
				strcat(inputname,DATAFILE_1);
        	        	strcat(outputname,OUTPUT_1);
			}
			if(file1ora==1) {
		        	strcat(inputname,DATAFILE_A);
                        	strcat(outputname,OUTPUT_A);
			}
			if(choice==2)
                                strcat(inputname,"HARD");
                        if(choice==3)
				strcat(inputname,"SIDE");
			if(choice==4)
				strcat(inputname,"SIDEHARD");
                        input=fopen(inputname,"r");

                        if(choice==2)
                                strcat(outputname,"HARD");
                        if(choice==3)
				strcat(outputname,"SIDE");
			if(choice==4)
				strcat(outputname,"SIDEHARD");

                        output=fopen(outputname,"w");

			/* now we take in the input file which has data that looks like:
			  sign, pe[0], pe[1], pe[2], pe[3], pg[0], pg[1], pg[2], pg[3]\n 
			*/
			while(1) {
				nscan1=fscanf(input,"%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf%c",
					&sign,&pe0,&pe1,&pe2,&pe3,&pg0,&pg1,&pg2,&pg3,&return1);
				if(nscan1==EOF)
					break;
				if(nscan1!=10) {
					printf("File error has occurred.\n");
					printf("nscan1==%d\n",nscan1);
					printf("Problem with %s\n",inputname);
					break;
				}
				pg0=pg0+Ge_res(uran,sptr);
				pe0=pe0+Ge_res(uran,sptr);
// changed Si resolution here				pe0=Si_res(pe0,uran,sptr);
				fprintf(output,"%lf, %lf, %d\n",pe0,pg0,sign);
			}
			fclose(input);
			fclose(output);
			output=fopen(outputname,"r");
			if(file1ora==0)
				histogram(output,K1count,average_1,filetype);
			if(file1ora==1)
				histogram(output,Kacount,average_a,filetype);
			fclose(output);
		}
	}
	histwrite(K1count,Kacount,average_1,average_a);
	return;
}

void histogram(FILE *file,int counts[NUM_BE_BINS][NUM_GE_BINS][NUM_GE],
	double average[NUM_BE_BINS][2][NUM_GE],int filetype)
{
	double BEdata,GEdata,BEtest,GEtest,dBE,dGE;
	int nscan1,i,j,sign;

	dBE=(BEHIST_END-BEHIST_START)/(double)NUM_BE_BINS;
	dGE=(GEHIST_END-GEHIST_START)/(double)NUM_GE_BINS;

	while(1) {
		nscan1=fscanf(file,"%lf, %lf, %d%*c",&BEdata,&GEdata,&sign);
		if(nscan1==EOF)
			break;
		if(nscan1!=3) {
                       	printf("Something went wrong in weighted_hist\n");
			break;
		}

		/* makes a 2D array with BE and GE histogram info */
		BEtest=BEHIST_START;
		for(i=0;i<NUM_BE_BINS;i++) {
			if((BEdata>=BEtest)&&(BEdata<BEtest+dBE)) {
/* start of new code */
				average[i][0][filetype]+=GEdata*sign;
				average[i][1][filetype]+=sign;
/* end of new code */
				GEtest=GEHIST_START;
				for(j=0;j<NUM_GE_BINS;j++) {
					if((GEdata>=GEtest)&&(GEdata<GEtest+dGE)) {
						counts[i][j][filetype]+=sign;
						i=NUM_BE_BINS+1; /* breaks from loop */
						j=NUM_GE_BINS+1;
					}
					GEtest+=dGE;
				}
			}
			BEtest+=dBE;
		}
	}
}

void histwrite(int K1[NUM_BE_BINS][NUM_GE_BINS][NUM_GE],int Ka[NUM_BE_BINS][NUM_GE_BINS][NUM_GE],
	double A1[NUM_BE_BINS][2][NUM_GE],double Aa[NUM_BE_BINS][2][NUM_GE])
{
	FILE *hist1,*hista,*ave1,*avea;
	int i,j,k;
	double dBE,dGE,BEbin,GEbin;

	hist1=fopen(HIST_1,"w");
	hista=fopen(HIST_A,"w");
	ave1=fopen(AVE_1,"w");
	avea=fopen(AVE_A,"w");
	dBE=(BEHIST_END-BEHIST_START)/(double)NUM_BE_BINS;
	dGE=(GEHIST_END-GEHIST_START)/(double)NUM_GE_BINS; 
	/* write this out as a text file that IGOR will read */
	fprintf(hist1,"IGOR\n");
	fprintf(hist1,"WAVES /D/O BE,GEbin,K1_1,K1_2,K1_3,K1_4\n");
	fprintf(hist1,"BEGIN\n");
	fprintf(hista,"IGOR\n");
	fprintf(hista,"WAVES /D/O BE,GEbin,Ka_1,Ka_2,Ka_3,Ka_4\n");
	fprintf(hista,"BEGIN\n");
	for(i=0;i<NUM_BE_BINS;i++) {
		for(j=0;j<NUM_GE_BINS;j++) {
			BEbin=BEHIST_START+((double)i+0.5)*dBE;
			GEbin=GEHIST_START+((double)j+0.5)*dGE;
			fprintf(hist1,"\t%.8lf\t%lf",1000.0*BEbin,GEbin);
			fprintf(hista,"\t%.8lf\t%lf",1000.0*BEbin,GEbin);
			for(k=0;k<NUM_GE;k++) {
				fprintf(hist1,"\t%d",K1[i][j][k]);
				fprintf(hista,"\t%d",Ka[i][j][k]);
			}
			fprintf(hist1,"\n");
			fprintf(hista,"\n");
		}
		fprintf(ave1,"%lf",1000.0*(BEbin-1.0)*mexp);
		fprintf(avea,"%lf",1000.0*(BEbin-1.0)*mexp);
		for(k=0;k<NUM_GE;k++) {
			if(A1[i][1][k]>0.0)
				fprintf(ave1,"\t%lf\t%lf\t%lf",
					A1[i][0][k]*mexp*1000.0,A1[i][1][k],(A1[i][0][k]*mexp/A1[i][1][k]-Egammaexp)*1000.0);
			else
				fprintf(ave1,"\t%lf\t%lf\t0.0",A1[i][0][k]*mexp*1000.0,A1[i][1][k]);
			if(Aa[i][1][k]>0.0)
				fprintf(avea,"\t%lf\t%lf\t%lf",
					Aa[i][0][k]*mexp*1000.0,Aa[i][1][k],(Aa[i][0][k]*mexp/A1[i][1][k]-Egammaexp)*1000.0);
			else
				fprintf(avea,"\t%lf\t%lf\t0.0",Aa[i][0][k]*mexp*1000.0,Aa[i][1][k]);
		}
		fprintf(ave1,"\n");
		fprintf(avea,"\n");
	}
	fprintf(hist1,"END\n");
	fprintf(hista,"END\n");
	fclose(hist1);
	fclose(hista);
	fclose(ave1);
	fclose(avea);
}

double Ge_res(double uran(long *), long *sptr)
{ 
        double res,sigma;
        /* For a Gaussian the relationship is sigma = FWHM/2.354 */
        sigma=GE_FWHM/2.354;
        res=sigma*GaussianDev(uran,sptr);
        return(res);
}

double Si_res(double E,double uran(long *), long *sptr)
{
	double res,response,randnum,Eprime,dE;
	int i;

	randnum=uran(sptr);
/*
	res=0.0;
	Eprime=0.0;
	dE=10.0/NUMCHANNELS;
	
	for(i=0;i<NUMCHANNELS;i++) {
		Eprime+=dE;
		response=(rc[0]*gaussian(E,Eprime)
			+ rc[1]*lowEtail(E,Eprime)
			+ rc[2]*highEplateau(E,Eprime)
			+ rc[3]*highEtail(E,Eprime))*dE;
		res+=response;
		printf("%lf\n",response);
	}
	printf("%lf	%lf\n",res,E);
*/
	response=0.0;
	Eprime=0.0;
	dE=10.0/NUMCHANNELS;
	
	for(i=0;i<NUMCHANNELS;i++) {
		Eprime+=dE;
		response+=(rc[0]*gaussian(E-1.0,Eprime)
			+ rc[1]*lowEtail(E-1.0,Eprime)
			+ rc[2]*highEplateau(E-1.0,Eprime)
			+ rc[3]*highEtail(E-1.0,Eprime))*dE;
		if(response>randnum)
			return(Eprime+1.0);
	}
	return(0.0);
}

double gaussian(double E,double Eprime)
{
	double gauss;
	gauss=1.0/(SIGMA*sqrt(2.0*PI))*exp(-(E-Eprime)*(E-Eprime)/(2.0*SIGMA*SIGMA));
	return(gauss);
}

double lowEtail(double E,double Eprime)
{
	double tail;
	tail=(1.0-erf((Eprime-E)/(SIGMA*sqrt(2.0))))/(2.0*E);
	return(tail);
}

double highEplateau(double E,double Eprime)
{
	double plateau;	
	plateau=(erf((Eprime-E)/(SIGMA*sqrt(2.0)))-erf((Eprime-E-COMPTON)/(SIGMA*sqrt(2.0))))/(2.0*COMPTON);
	return(plateau);
}

double highEtail(double E,double Eprime)
{
	double hightail;	
	hightail=((Eprime-E)*(erf((Eprime-E)/(SIGMA*sqrt(2.0)))-2.0*erf((Eprime-E-COMPTON)/(SIGMA*sqrt(2.0)))
			+ erf((Eprime-E-2.0*COMPTON)/(SIGMA*sqrt(2.0))))
		+ 2.0*COMPTON*(erf((Eprime-E-COMPTON)/(SIGMA*sqrt(2.0)))-erf((Eprime-E-2.0*COMPTON)/(SIGMA*sqrt(2.0)))) 
		+ 2.0*SIGMA*SIGMA*(gaussian(E,Eprime)-2.0*gaussian(E+COMPTON,Eprime)
			+ gaussian(E+2.0*COMPTON,Eprime)))/(2.0*COMPTON*COMPTON);
	return(hightail);
}
