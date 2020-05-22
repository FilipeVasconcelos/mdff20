#ifndef CONFIG_H
#define CONFIG_H
#include "constants.h"
#include "mpi_mdff.h"
/* ****************************************************************************/
/*                              global parameters                             */
/* ****************************************************************************/
int                    nion; 
double              onenion; /* 1/nion     */
double         onethirdnion; /* 1/(3*nion) */
int                   ntype;
double    massit [NTYPEMAX];
int       nionit [NTYPEMAX];
char*     atypit [NTYPEMAX];

char  configname[MAX_LEN+1];

double          *rx,*ry,*rz;  
double          *vx,*vy,*vz;  
double          *fx,*fy,*fz;  
double       *rxs,*rys,*rzs;  
double              *massia;  
double          *invemassia;  
char               **atypia;  
int                  *typia; 
double       tau_nonb[3][3];
double                 rhoN; /*density number:   rhoN=nion/V  (Å^-3)*/ 

/* ****************************************************************************/
/*                                prototypes                                  */
/* ****************************************************************************/
void          init_config();
void         alloc_config();
void          free_config();
void          info_config();
void sample_config(int key);
int          write_config();
int           read_config();
#endif
