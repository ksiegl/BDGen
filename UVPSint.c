#include "UVPSpline.h"

#include "UVPSint.h"

#include "UVPSloc.h"



double UVPSint(double x0,double x1,const UVPSpline *sptr)

{

	const UVPSData *d0ptr,*d1ptr,*dptr0,*dptr1;

	double fval,u0,u1,w0,w1;

	long p;

	const double *aptr;

	

	fval=0.0;

	

	if (sptr) {

		if (x0<x1) {

			if ((dptr0=UVPSloc(x0,sptr->dptr0,sptr->dptr1))&&(dptr1=UVPSloc(x1,sptr->dptr0,sptr->dptr1))) {

				for (d0ptr=dptr0;d0ptr<=dptr1;++d0ptr) {

					d1ptr=d0ptr+1;

					u0=(x0>d0ptr->x)?x0-d0ptr->x:0.0;

					if (d1ptr<sptr->dptr1) {

						u1=(x1<d1ptr->x)?x1-d0ptr->x:d1ptr->x-d0ptr->x;

					} else {

						u1=x1-d0ptr->x;

					}

					p=0;

					w0=u0;

					w1=u1;

					for (aptr=d0ptr->aptr0;aptr<d0ptr->aptr1;++aptr) {

						++p;

						fval+=(*aptr)*(w1-w0)/((double)p);

						w0*=u0;

						w1*=u1;

					}

				}

			}

		}

	}

	

	return fval;

}

