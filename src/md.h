#ifndef MD_H
#define MD_H
/* ****************************************************************************/
/*                              global parameters                             */
/* ****************************************************************************/
int                   istep; /* current step */
int                    npas; /* number of md steps */
int                  nprint; /* stdout printing period */
int                  fprint; /* osziff printing period */
int                  cprint; /* contff printing period */
int                  tprint; /* trajff printing period */
int                  tstart; /* trajff starting step   */
int                  nequil; /* equilibration interval (rescaling) */
int                 nequilT; /* equilibration period (rescaling) */
double                 temp; /* temperature */
double                   dt; /* time step */
double        tauTberendsen; /* berendsen thermostat period */
int                 egrator; /* int version of */
/*****************************************************************************/
/*                  Nose-Hoover Chain related parameters                     */
/*****************************************************************************/
double       timesca_thermo; /* time scale of thermostat */
double         timesca_baro; /* time scale of barostat */
int          nhc_yosh_order; /* order of the yoshida integrator */
int               nhc_mults; /* number of multiple timesteps  */
int                   nhc_n; /* length of the Nose-Hoover chain */
double           *vxi,  *xi; /* thermostat coordinates coupled to the particules */
double          *vxib, *xib; /* barostat coordinates coupled to the volume */

#define ALLWD_YOSH_PARAM 5
int yosh_allowed[ALLWD_YOSH_PARAM];

#define ALLWD_INTEGRATOR_STR 3
char allwd_integrator[ALLWD_INTEGRATOR_STR][MAX_LEN+1];
#define ALLWD_RESCALE_INTEGRATOR 1
int allwd_rescale_integrator[ALLWD_RESCALE_INTEGRATOR];


#include <stdbool.h>
bool             lleapequi; /* true if integrator == nve-lf and nequil > 0 */
bool       rescale_allowed;
/* ****************************************************************************/
/*                                prototypes                                  */
/* ****************************************************************************/
int read_md(char* controlfn);
void init_md(char* controlfn);
void info_md();
void run_md();
void alloc_md();
void free_md();
int check_md();
#endif
