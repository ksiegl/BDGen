#include "UVPSpline.h"

#include "UVPSval.h"

#include "UVPSloc.h"



double UVPSval(double x,const UVPSpline *sptr)

{

	const UVPSData *dptr;

	double fval,u,w;

	const double *aptr;

	

	fval=0.0;

	

	if (sptr) {

		if (dptr=UVPSloc(x,sptr->dptr0,sptr->dptr1)) {

			u=x-dptr->x;

			w=1.0;

			for (aptr=dptr->aptr0;aptr<dptr->aptr1;++aptr) {

				fval+=w*(*aptr);

				w*=u;

			}

		}

	}

	

	return fval;

}

