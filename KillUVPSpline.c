#include <stdlib.h>

#include "UVPSpline.h"

#include "KillUVPSpline.h"



UVPSpline *KillUVPSpline(UVPSpline *sptr)

{

	UVPSData *dptr;

	

	if (sptr) {

		if (sptr->dptr0) {

			for (dptr=sptr->dptr0;dptr<sptr->dptr1;++dptr) {

				if (dptr->aptr0) {

					free(dptr->aptr0);

					dptr->aptr0=NULL;

					dptr->aptr1=NULL;

				}

			}

			free(sptr->dptr0);

			sptr->dptr0=NULL;

			sptr->dptr1=NULL;

		}

		free(sptr);

	}

	

	return NULL;

}

