#ifndef CONSTANTS_H
#define CONSTANTS_H

#define MAX_LEN 80
#define NTYPEMAX 16
#define NPAIRSMAX 136 
#define LOWERCASE 97

#define PI        3.1415926535897932385  /* pi  */
#define TPI       6.2831853071795864770  /* 2pi */
#define xx " "
#define xx3 "   "
#define xx6 xx3 xx3
#define xx9 xx6 xx xx
#define xx11 xx9 xx3
#define ff "%15.8f "
#define ee "%15.8e "
#define FF "%19.12f "
#define EE "%19.12e "
#define EE2 "%19.12e %19.12e"
#define ff3 "%15.8f %15.8f %15.8f "
#define ee3 "%15.8e %15.8e %15.8e "
#define FF3 "%19.12f %19.12f %19.12f "
#define EE3 "%19.12e %19.12e %19.12e "
#define ff5 "%15.8f %15.8f %15.8f %15.8f %15.8f "
#define ee6 ee3 ee3
#define ee9 ee6 ee3

#define tt3 "%14.7e %14.7e %14.7e "
/*
--------------------------------------------------
  Main constants and units definition of the code
--------------------------------------------------

Units of MDFF in input/output :
length        : angstom
energy        : eV
mass          : atomic mass unit
time          : ps
temperature   : Kelvin
pressure      : GPa
charge        : atomic charge ( proton charge )
dipole moment : Debye
efg           : V  / angstom**2
velocity      : ( atomicmassunit / eV ) ** 0.5
forces        : eV / angtrom
density       : g/cm^3

Internal units:
length        : angstom
energy        : eV
mass          : atomic mass unit
time          : angstrom * ( atomicmassunit / eV ) ** 0.5
temperature   : eV
pressure      : ( eV / angstrom**3 )
charge        : atomic charge ( proton charge )
dipole moment : angstrom * atomiccharge
efg           : atomiccharge / angstrom ** 3
velocity      : ( atomicmassunit / eV ) ** 0.5
forces        : eV / angtrom
density       : atomic mass unit / angstrom^3
1/4piepsilon0 : eV  * angstrom / atomiccharge **2 == 1.0
*/

char allwd_FT_str[16][MAX_LEN+1];

double radian  ; /* 180/pi */
double onethird; /* 1/3    */
double piroot  ;

double boltz_unit;  /* boltzmann constant ( energy in eV) */
double time_unit ;  /* unit of time picosecond => angstrom * ( atomicmassunit / eV ) ** 0.5 */
double press_unit;  /* GPa                     => internal unit of pressure ( eV / angstrom**3) */
double rho_unit  ;  /* g/cm^3                  => internal unit of density (atomicmassunit/angstrom^3)*/
double coul_unit ;  /*  1 / 4pi epsilon0  in eV */


/* ****************************************************************************/
/*                                prototypes                                  */
/* ****************************************************************************/
void gen_constants();
void reduced_units();
#endif
