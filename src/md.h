#ifndef MD_H
#define MD_H
/* ****************************************************************************/
/*                              global parameters                             */
/* ****************************************************************************/
int                   istep; /* current step */
int                    npas; /* number of md steps */
int                  nprint; /* stdout printing period */
int                  fprint; /* osziff printing period */
int                  nequil; /* equilibration period (rescaling) */
double                 temp; /* temperature */
double                   dt; /* time step */  
double        tauTberendsen; /* Berendsen thermostat period */
int                 egrator; /* int version of */

#define ALLWD_INTEGRATOR_STR 2
char allwd_integrator_str[ALLWD_INTEGRATOR_STR][MAX_LEN+1];

#include <stdbool.h>
bool             lleapequi; /* true if integrator == nve-lf and nequil > 0 */

/* ****************************************************************************/
/*                                prototypes                                  */
/* ****************************************************************************/
int read_md(char* controlfn);
void init_md(char* controlfn);
void info_md();
void run_md();
#endif
