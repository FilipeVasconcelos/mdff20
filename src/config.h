#ifndef CONFIG_H
#define CONFIG_H
#include "constants.h"
#include "mpi_mdff.h"
int nion; double onenion;
int ntype;

char configname[MAX_LEN+1];

double *rx,*ry,*rz;  
double *vx,*vy,*vz;  
double *fx,*fy,*fz;  
double *rxs,*rys,*rzs;  
double *massia;  
double *invemassia;  
//char   *atype[MAX_LEN+1];  
char   **atype;  
int     *itype; 
double tau_nonb[3][3];


int read_config();
void init_config();
void alloc_config();
#endif
