#include <stdlib.h>

#include "UVPSpline.h"

#include "UVPSloc.h"



const UVPSData *UVPSloc(double x,const UVPSData *dptr0,const UVPSData *dptr1)

{

	const UVPSData *dptr;

	long n;

	

	if (dptr0<dptr1) {

		if (x<dptr0->x) {

			dptr=dptr0;

		} else if (x>=(dptr1-1)->x) {

			dptr=dptr1-1;

		} else {

			while ((n=(long)(dptr1-dptr0))>1) {

				dptr=dptr0+n/2;

				if (x>=dptr->x) {

					dptr0=dptr;

				} else {

					dptr1=dptr;

				}

			}

			dptr=dptr0;

		}

	} else {

		dptr=NULL;

	}

	

	return dptr;

}

