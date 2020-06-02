#ifndef THERMO_H
#define THERMO_H
double        e_kin; /* kinetic energy */
double        e_nvt; /* Nose Hoover thermostat "energy" (conserved quantity) */
double         u_lj; /* nmlj potential energy */
double       u_coul; /* coulombic potential energy */
double      pvir_lj; /* Pvirial nmlj */ 
double    pvir_coul; /* Pvirial coulombic */
double tau_lj[3][3]; /* nmlj stress tensor */


/* ****************************************************************************/
/*                                prototypes                                  */
/* ****************************************************************************/
void info_thermo();
#endif
