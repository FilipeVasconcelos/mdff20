
double     *efia;  /* electric field at ion position */
double    *efgia;  /* electric field gradient at ion position */


double      *qia;  /* charge     at ion */
double     *muia;  /* dipole     at ion */
double  *thetaia;  /* quadrupole at ion */ 


/* ****************************************************************************/
/*                                prototypes                                  */
/* ****************************************************************************/
void set_autoES();

void multipole_ES(double *q, double (*mu)[3], double (*theta)[3][3],double *u);

void multipole_ES_dir(double *q, double (*mu)[3], double (*theta)[3][3],
                      double *u_dir  , double (*ef_dir)[3], double (*efg_dir)[3][3], 
                      double *fx_dir, double *fy_dir, double *fz_dir , double tau_dir[3][3]);
void multipole_ES_rec(double *q, double (*mu)[3], double (*theta)[3][3],
                      double *u_rec  , double (*ef_rec)[3], double (*efg_rec)[3][3],
                      double *fx_rec, double *fy_rec, double *fz_rec , double tau_rec[3][3]);


