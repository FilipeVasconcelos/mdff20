#ifndef KSPACE_H
#define KSPACE_H

#include "constants.h"
#include "mpi_mdff.h"
struct KMESH kcoul;

typedef struct KMESH {
int                      nk; /* total number of kpoints                                          */
int                 kmax[3]; /* nb of kpts in each recip. direction in ewald sum.                */
double          *kx,*ky,*kz; /* kpoints mesh  array dimension [nk]                               */
double                  *kk; /* k module                                                         */
double                  *Ak; /* Ak quantity in ewald  exp( - kk * 0.25 / alpha2 ) / kk           */
double                *kcoe; /* kcoe quantity in ewald   2.0 * ( 1.0 / kk + 1.0 / alpha2 / 4.0 ) */
char   meshlabel[MAX_LEN+1]; /* giving a name to kmesh                                           */
DEC                  kptDec; /* k-point decomposition                                            */
double              **ckria; /* cos( k r )                                                       */
double              **skria; /* sin( k r )                                                       */
double              *rhon_R; /* Re ( charge density )                                            */
double              *rhon_I; /* Im ( charge density )                                            */
double                 *str; /* (rhonk_R*rhonk_R + rhonk_I*rhonk_I) * Ak                         */
} KMESH;

void init_kspace();
void free_kspace();
void set_param_kmesh(KMESH *km,double alpha);
void reorder_kmesh(KMESH *km);
void struct_fact_rhon(double *q, double (*mu)[3], double (*theta)[3][3],bool lqchtask, bool ldiptask);
#endif
