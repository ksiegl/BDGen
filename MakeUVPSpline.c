#include <stdlib.h>
#include "UVPSpline.h"
#include "MakeUVPSpline.h"
#include "KillUVPSpline.h"

UVPSpline *MakeUVPSpline(long len)
{
	UVPSpline *sptr;
	UVPSData *dptr;
	
	if (len<1) return NULL;
	
	if (!(sptr=(UVPSpline *)malloc(sizeof(UVPSpline)))) 
		return NULL;
	sptr->dptr0=NULL;
	sptr->dptr1=NULL;
	
	if (!(sptr->dptr0=(UVPSData *)malloc(len*sizeof(UVPSData)))) 
		return KillUVPSpline(sptr);
	sptr->dptr1=sptr->dptr0+len;
	
	for (dptr=sptr->dptr0;dptr<sptr->dptr1;++dptr) {
		dptr->aptr0=NULL;
		dptr->aptr1=NULL;
	}
	return sptr;
}
