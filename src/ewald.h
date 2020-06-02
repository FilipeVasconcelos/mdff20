
double     *efia;  /* electric field at ion position */
double    *efgia;  /* electric field gradient at ion position */


double      *qia;  /* charge     at ion */
double     *muia;  /* dipole     at ion */
double  *thetaia;  /* quadrupole at ion */ 

/* ****************************************************************************/
/*                                prototypes                                  */
/* ****************************************************************************/
void multipole_ewald_sum();
void multipole_ewald_sum_direct(double u_dir  , double (*ef_dir)[3], double (*efg_dir)[3][3], 
                                double *fx_dir, double *fy_dir, double *fz_dir , double tau_dir[3][3]);
void multipole_ewald_sum_reciprocal(double u_rec  , double (*ef_rec)[3], double (*efg_rec)[3][3],
                                    double *fx_rec, double *fy_rec, double *fz_rec , double tau_rec[3][3]);


