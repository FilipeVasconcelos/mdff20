#ifndef PIM_H
#define PIM_H

#define MAX_EXTRAPOLATE_ORDER 5
int extrapolate_order;
int min_scf_pim_iter;
int max_scf_pim_iter;
double omegakO;

#include <stdbool.h>
double conv_tol_ind;

#define ALLWD_ALGO_PIM_STR 2
char allwd_algo_pim[ALLWD_ALGO_PIM_STR][MAX_LEN+1];
int algo_pim;
#define ALLWD_ALGO_EXTRAPOLATE_DIPOLE_STR 2
char allwd_algo_extrapolate_dipole[ALLWD_ALGO_EXTRAPOLATE_DIPOLE_STR][MAX_LEN+1];
int algo_extrapolate_dipole;


/* ****************************************************************************/
/*                                prototypes                                  */
/* ****************************************************************************/
void info_pim();
void init_pim(char * controlfn);
void momentpolaSCF(double (*mu)[3],double *upol);
void momentpolaSCFkO(double (*mu_ind)[3],double *upol);
void induced_moment(double (*mu_ind)[3], double (*ef)[3]);
void induced_moment_inner_kO(double (*mu_ind)[3], double (*ef_ext)[3]);
void extrapolate_dipole_aspc(double (*mu_ind)[3] , double (*ef)[3], int key );
double get_rmsd_mu(double (*mu)[3], double (*ef)[3]);
void find_min(double (*ef0)[3], double (*ef1)[3],double (*mu0)[3], double (*mu1)[3]) ;
#endif
