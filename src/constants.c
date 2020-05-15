#include<stdlib.h>
#include<stdio.h>
#include "constants.h"

void gen_constants()
{

    boltz_unit = 8.61734229648141E-05F;  // boltzmann constant ( energy in eV ) 
    time_unit  = 98.2269514139276F    ;  // unit of time picosecond => angstrom * ( atomicmassunit / eV ) ** 0.5
    press_unit = 0.00624150964712042F ;  // GPa => internal unit of pressure (eV/A^3)
    rho_unit   = 1.660538782F         ;

    radian   = 180. / PI;     
    onethird = 1./3.0;

    for (int it=0;it<NTYPEMAX;it++)
    {
        atypei[it] = malloc((MAX_LEN+1)*sizeof (char*));
        mass[it]=1.0;
        natmi[it]=0;
    }

}


void reduced_units()
{
    boltz_unit = 1.0;
    time_unit  = 1.0;
    press_unit = 1.0;
    rho_unit   = 1.0;
}
