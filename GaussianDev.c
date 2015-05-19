#include <math.h>
#include "GaussianDev.h"

static int iset=0;
static double gset;

/* This routine comes from "Numerical Recipes in C" p. 289-290 */
/* Returns a normally distributed deviate with zero mean and unit variance */

double GaussianDev(double uran(long *),long *sptr)
{
	float fac,rsq,v1,v2;

	/* since we don't have an extra deviate handy */
	if(iset==0) {
		do {
			/* pick two uniform numbers in the square extending 
			   from -1 to +1 in each direction and see if they
			   are in the unit circle */
			v1=2.0*uran(sptr)-1.0;
			v2=2.0*uran(sptr)-1.0;
			rsq=v1*v1+v2*v2;
		} while(rsq>=1.0 || rsq==0.0); /* and if they are not
						  try again */
		fac=sqrt(-2.0*log(rsq)/rsq);
		/* Now make the Box-Muller transformation to get two normal
		   deviates.  Return one and save the other for next time. */
		gset=v1*fac;
		iset=1;				/* set flag */
		return(v2*fac);
	} else {	
		/* we have an extra deviate handy, so unset the flag and
		   return it */
		iset=0;
		return(gset);
	}
}
