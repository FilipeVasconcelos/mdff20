#ifndef THERMO_H
#define THERMO_H
double            e_kin; /* kinetic energy */
double            e_nvt; /* Nose Hoover thermostat "energy" (conserved quantity) */
double             u_lj; /* nmlj potential energy */
double         u_bhmftd; /* BHMFTD potential energy */
double           u_coul; /* coulombic potential energy */
double            u_pol; /* polarizabily energy */
double          pvir_lj; /* Pvirial nmlj */
double      pvir_bhmftd; /* Pvirial BHMFTD */
double        pvir_coul; /* Pvirial coulombic */
double     tau_lj[3][3]; /* nmlj stress tensor */
double tau_bhmftd[3][3]; /* BHMFTD stress tensor */
double   tau_coul[3][3]; /* electrostatic stress tensor */


/* ****************************************************************************/
/*                                prototypes                                  */
/* ****************************************************************************/
void info_thermo();
#endif
