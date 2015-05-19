#include <stdlib.h>
#include "UVPSpline.h"
#include "CompUVPASpline.h"
#include "MakeUVPSpline.h"
#include "KillUVPSpline.h"

UVPSpline *CompUVPASpline(long len,const double *x,const double *y)
{
	UVPSpline *sptr;
	UVPSData *dptr;
	long i,i0,i1,i2,i3,j,k;
	double x0,x1,x2,x3,y0,y1,y2,y3,dx,dy,dy0,dy1,dy2,dy3,yp,yp0,yp1;
	double d,sx,sy,sxx,sxy,b0,b1,vol,pe,wt,eps;
	double smpef,smwtf,smpei,smwti;
	
	if (len<5) return NULL;
	
	for (i=1;i<len;++i) {
		if (x[i]<=x[i-1]) return NULL;
	}
	
	if (!(sptr=MakeUVPSpline(len))) return NULL;
	
	dptr=sptr->dptr0;
	
	for (i=0;i<len-1;++i) {
		if (!(dptr->aptr0=(double *)malloc(4*sizeof(double)))) 
			return KillUVPSpline(sptr);
		dptr->aptr1=dptr->aptr0+4;
		for (j=0;j<2;++j) {
			i0=i+j;
			x0=x[i0];
			y0=y[i0];
			smpef=0.0;
			smwtf=0.0;
			smpei=0.0;
			smwti=0.0;
			for (k=0;k<4;++k) {
				if (k==0) {
					i1=i0-3;
					i2=i0-2;
					i3=i0-1;
				} else if (k==1) {
					i1=i0+1;
				} else if (k==2) {
					i2=i0+2;
				} else {
					i3=i0+3;
				}
				if ((i1>=0)&&(i1<len)&&(i2>=0)&&(i2<len)&&(i3>=0)&&(i3<len)) {
					x1=x[i1]-x0;
					x2=x[i2]-x0;
					x3=x[i3]-x0;
					y1=y[i1]-y0;
					y2=y[i2]-y0;
					y3=y[i3]-y0;
					d=x1*x2*x3*(x2-x1)*(x3-x2)*(x3-x1);
					pe=(x2*x2*x3*x3*(x3-x2)*y1+x3*x3*x1*x1*(x1-x3)*y2+x2*x2*x1*x1*(x2-x1)*y3)/d;
					sx=x1+x2+x3;
					sy=y1+y2+y3;
					sxx=x1*x1+x2*x2+x3*x3;
					sxy=x1*y1+x2*y2+x3*y3;
					d=4.0*sxx-sx*sx;
					b0=(sxx*sy-sx*sxy)/d;
					b1=(4.0*sxy-sx*sy)/d;
					dy0=(-b0);
					dy1=y1-b0-b1*x1;
					dy2=y2-b0-b1*x2;
					dy3=y3-b0-b1*x3;
					vol=dy0*dy0+dy1*dy1+dy2*dy2+dy3*dy3;
					eps=1.0e-12*(y[i0]*y[i0]+y[i1]*y[i1]+y[i2]*y[i2]+y[i3]*y[i3]);
					if (vol>eps) {
						wt=1.0/(vol*sxx);
						smpef+=pe*wt;
						smwtf+=wt;
					} else {
						smpei+=pe;
						smwti+=1.0;
					}
				}
			}
			if (smwti<0.5) {
				yp=smpef/smwtf;
			} else {
				yp=smpei/smwti;
			}
			if (j==0) {
				yp0=yp;
			} else {
				yp1=yp;
			}
		}
		dx=x[i+1]-x[i];
		dy=y[i+1]-y[i];
		dptr->x=x[i];
		dptr->aptr0[0]=y[i];
		dptr->aptr0[1]=yp0;
		yp1-=yp0;
		yp0-=dy/dx;
		dptr->aptr0[2]=(-(3.0*yp0+yp1)/dx);
		dptr->aptr0[3]=(2.0*yp0+yp1)/(dx*dx);
		++dptr;
	}
	
	if (!(dptr->aptr0=(double *)malloc(2*sizeof(double)))) return
KillUVPSpline(sptr);
	dptr->aptr1=dptr->aptr0+2;
	
	x0=x[len-1];
	x1=x[len-2]-x0;
	x2=x[len-3]-x0;
	x3=x[len-4]-x0;
	y0=y[len-1];
	y1=y[len-2]-y0;
	y2=y[len-3]-y0;
	y3=y[len-4]-y0;
	d=x1*x2*x3*(x2-x1)*(x3-x2)*(x3-x1);
	
	dptr->x=x0;
	dptr->aptr0[0]=y0;
	dptr->aptr0[1]=(x2*x2*x3*x3*(x3-x2)*y1+x3*x3*x1*x1*(x1-x3)*y2+x2*x2*x1*x1*(x2-x1)*y3)/d;
	
	return sptr;
}
