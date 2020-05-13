#ifndef CONFIG_H
#define CONFIG_H
#include "constants.h"
int nion; double onenion;
int ntype;

double *rx,*ry,*rz;  
double *vx,*vy,*vz;  
double *fx,*fy,*fz;  
double *rxs,*rys,*rzs;  
double *massia;  
double *invemassia;  
char   **atype;  
int     *itype; 
double tau_nonb[3][3];

int read_config();
void init_config();
#endif
