#include "BD2.h"
#include <math.h>
#include "c3ranradius.h"

void c3ranradius(double *cxptr,double *cyptr,double *czptr,double uran(long *),long *sptr)
{
	double rsq;
	
	for (;;) {
		*cxptr=2.0*uran(sptr)-1.0;
		*cyptr=2.0*uran(sptr)-1.0;
		*czptr=2.0*uran(sptr)-1.0;
		rsq=(*cxptr)*(*cxptr)+(*cyptr)*(*cyptr)+(*czptr)*(*czptr);
		if ((rsq>0.0000000)&&(rsq<1.00000000)) {
			*cxptr*=TRAP_FWHM;
			*cyptr*=TRAP_FWHM;
			*czptr*=TRAP_FWHM;
			return;
		}
	}
}
