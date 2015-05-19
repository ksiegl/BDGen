#include <float.h>
#include "ran2.h"

#define IM1 2147483563
#define IM2 2147483399
#define AM1 (1.0/(1.0*IM1))
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNOMAX 1.0-DBL_EPSILON

static long seed2=123456789;
static long seed0=0;
static long data[NTAB];

double ran2(long *seed1)
{
	int i;
	long k;
	double rno;
	
	if (((*seed1)<=0)||(seed0==0)) {
		if ((*seed1)<0) (*seed1)=-(*seed1);
		else if ((*seed1)==0) (*seed1)=1;
		seed2=(*seed1);
		for (i=NTAB+7;i>=0;--i) {
			k=(*seed1)/IQ1;
			(*seed1)=IA1*((*seed1)-k*IQ1)-k*IR1;
			if ((*seed1)<0) (*seed1)+=IM1;
			if (i<NTAB) data[i]=(*seed1);
		}
		seed0=data[0];
	}
	
	k=(*seed1)/IQ1;
	(*seed1)=IA1*((*seed1)-k*IQ1)-k*IR1;
	if ((*seed1)<0) (*seed1)+=IM1;
	k=seed2/IQ2;
	seed2=IA2*(seed2-k*IQ2)-k*IR2;
	if (seed2<0) seed2+=IM2;
	i=seed0/NDIV;
	seed0=data[i]-seed2;
	data[i]=(*seed1);
	if (seed0<1) seed0+=IMM1;
	
	return ((rno=AM1*seed0)<RNOMAX)?rno:RNOMAX;
}
