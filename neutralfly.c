#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "BD2.h"

#define LIGHTSPEED 2.997924562e+5 /* c in mm/us */
int main(void)
{
	FILE *input,*output;	
	int i,nscan1;
	double locx,locy,locz,KE_SimIon;
	double az,el,norm,velx,vely,velz,time,radiusy,radiusz,radius;
	char inname[25],outname[25];

        printf("Input data file ==>");
        scanf("%s",inname);
        printf("Output data file ==>");
        scanf("%s",outname);

        input=fopen(inname,"r");
        output=fopen(outname,"w");

	i=0;
	fprintf(output,"------ Begin Next Fly'm ------\n");
        while(1) {
		i++;
                nscan1=fscanf(input,"%*d, %*lf, %*d, %lf, %lf, %lf, %lf, %lf, %lf, %*d, %*d%*c",
			&locx,&locy,&locz,&az,&el,&KE_SimIon);
                if(nscan1==EOF)
                        break;
                if(nscan1!=6) {
                        printf("Something went wrong in reading the file.\n");
                        break;
		}
		velx=sqrt(2.0*KE_SimIon/(M1exp*uexp))*cos(az*PI/180.0)*cos(el*PI/180.0)*LIGHTSPEED;
		vely=sqrt(2.0*KE_SimIon/(M1exp*uexp))*sin(az*PI/180.0)*cos(el*PI/180.0)*LIGHTSPEED;
		velz=sqrt(2.0*KE_SimIon/(M1exp*uexp))*sin(el*PI/180.0)*LIGHTSPEED;

		time=(MCP_X-locx)/velx;
		radiusy=locy-CENTER_Y+vely*time;
		radiusz=locz-CENTER_Z+velz*time;
		radius=sqrt(radiusy*radiusy+radiusz*radiusz);
		radiusy+=CENTER_Y;
		radiusz+=CENTER_Z;
		if(radius<=MCP_R && time>0.0)
			fprintf(output,"%d,%lf,%lf,%lf,%lf,%lf\n",i,time,MCP_X,radiusy,radiusz,KE_SimIon);
		else
			fprintf(output,"%d,%lf,%lf,%lf,%lf,%lf\n",i,time,MCP_X+1000,radiusy,radiusz,KE_SimIon);	
	}              
}
