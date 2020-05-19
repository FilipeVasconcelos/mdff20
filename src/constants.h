#ifndef CONSTANTS_H
#define CONSTANTS_H
#define NTYPEMAX 16
#define MAX_LEN 80
#define LOWERCASE 97
#define PI        3.1415926535897932385  /* pi  */
#define TPI       6.2831853071795864770  /* 2pi */ 

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
forces        : eV / angtrom
density       :  atomic mass unit / angstrom^3
1/4piepsilon0 : eV  * angstrom / atomiccharge **2 == 1.0
*/

char allowed_Tstring[8][MAX_LEN+1];
char allowed_Fstring[8][MAX_LEN+1];


double radian  ; /* 180/pi */
double onethird; /* 1/3    */

double boltz_unit;  /* boltzmann constant ( energy in eV) */
double time_unit ;  /* unit of time picosecond => angstrom * ( atomicmassunit / eV ) ** 0.5 */
double press_unit;  /* GPa                     => internal unit of pressure ( eV / angstrom**3) */
double rho_unit  ;  /* g/cm^3                  => internal unit of density (atomicmassunit/angstrom^3)*/

double mass   [NTYPEMAX];
int    natmi  [NTYPEMAX];
char*  atypei [NTYPEMAX];

void gen_constants();
#endif
