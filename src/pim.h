#ifndef PIM_H
#define PIM_H
#include <stdbool.h>
double conv_tol_ind;

#define ALLWD_ALGO_PIM_STR 2
char allwd_algo_pim[ALLWD_ALGO_PIM_STR][MAX_LEN+1];
int algo_pim;
#define ALLWD_ALGO_EXTRAPOLATE_DIPOLE_STR 2
char allwd_algo_extrapolate_dipole[ALLWD_ALGO_EXTRAPOLATE_DIPOLE_STR][MAX_LEN+1];
int algo_extrapolate_dipole;

int extrapolate_order;

void info_pim();
void momentpolaSCF(double (*mu)[3]);
#endif
