
void BDA021Z11(double *w1,int *w_sign,double *pe,double *pv,double *p1,double *pg,double pgarray[][4],double pcearray[][4],double *pn,
	       double *pe_exp,double *pv_exp,double *p1_exp,double *pg_exp,double *pn_exp,
	       double *loc,double *loc_n,
	       double *beta_rad,double *MCP_hit,double uran(long *),long *sptr,double M0exp,double M1exp,double M2exp,double Sn_exp,double abv,
		   int numGammas,double EGammas, double *Gammas,int cenum, double *ce, double Neutron1, int choice, double maxdist);
double dGamma1(double E,double cos,double E0,double M,double cpex,double cpey,double cpez,double cpvx,
	       double cpvy,double cpvz,int *w_sign);