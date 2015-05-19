#include <math.h>
#include "c3ran.h"

void c3ran(double *cxptr,double *cyptr,double *czptr,double uran(long *),long *sptr)
{
	double r,rsq;
	
	for (;;) {
		*cxptr=2.0*uran(sptr)-1.0;
		*cyptr=2.0*uran(sptr)-1.0;
		*czptr=2.0*uran(sptr)-1.0;
		rsq=(*cxptr)*(*cxptr)+(*cyptr)*(*cyptr)+(*czptr)*(*czptr);
		if ((rsq>0.0000000)&&(rsq<1.00000000)) {
			r=sqrt(rsq);
			*cxptr/=r;
			*cyptr/=r;
			*czptr/=r;
			return;
		}
	}
}
