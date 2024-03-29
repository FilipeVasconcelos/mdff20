#ifndef CONFIG_H
#define CONFIG_H
#include "constants.h"
#include "mpi_mdff.h"
#include "pim.h"
/* ****************************************************************************/
/*                              global parameters                             */
/* ****************************************************************************/
int                      nion;
double                onenion; /* 1/nion     */
double           onethirdnion; /* 1/(3*nion) */
int                     ntype;
char    configname[MAX_LEN+1];
double                   rhoN; /*density number:   rhoN=nion/V  (Å^-3)*/
/* ****************************************************************************/
/*                            on type quantities                              */
/* ****************************************************************************/
double       massit[NTYPEMAX];
int        nionit[NTYPEMAX+1]; /* last one is nion */
double invenionit[NTYPEMAX+1]; /* 1/nionit */
char*        atypit[NTYPEMAX]; /* string type */
double          qit[NTYPEMAX]; /* charge type */
double     dipit[NTYPEMAX][3]; /* dipole type */
double quadit[NTYPEMAX][3][3]; /* quadrupole type */
double  polit[NTYPEMAX][3][3]; /* dipole polarizability tensor if lpolar[it] = true*/
/* dipole damping */
bool   lpoldamping[NTYPEMAX][NTYPEMAX][NTYPEMAX];
/*
  lpoldamping[i][j][k] for type i and interaction bewteen type j and k
*/
double  pol_damp_b[NTYPEMAX][NTYPEMAX][NTYPEMAX];
double  pol_damp_c[NTYPEMAX][NTYPEMAX][NTYPEMAX];
int     pol_damp_k[NTYPEMAX][NTYPEMAX][NTYPEMAX];

/* ****************************************************************************/
/*                           on ion quantities                                */
/* ****************************************************************************/
double                            *rx,*ry,*rz; /* postion of ion */
double                            *vx,*vy,*vz; /* velocitie of ion */
double                            *fx,*fy,*fz; /* forces on ion */
double                         *rxs,*rys,*rzs; /* save position */
double                                *massia; /* mass of ion */
double                            *invemassia; /* inverse of masse of ion */
char                                 **atypia; /* string type of ion */
int                                    *typia; /* type of ion */
double                                   *qia; /* charge of ion */
double                            (*dipia)[3]; /* dipole of ion */
double (*dipia_ind)[MAX_EXTRAPOLATE_ORDER+1][3]; /* save induced dipole of ion */
double                        (*quadia)[3][3]; /* quadrupole of ion */
double                         (*polia)[3][3]; /* polarizability tensor on ion */
double                     (*invepolia)[3][3]; /* inverse polarizability tensor on ion */
/* ****************************************************************************/
/*                                prototypes                                  */
/* ****************************************************************************/
void          init_config();
void         alloc_config();
void          free_config();
void          info_config();
void sample_config(int key);
void sample_(char *label,double *ax, double *ay, double *az);
int          write_config(char* filename, char opt);
int           read_config();
void          read_next_traj(FILE *fp);
#endif
