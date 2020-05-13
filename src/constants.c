#include<stdlib.h>
#include<stdio.h>
#include "constants.h"

void gen_constants()
{
    TPI=2.0*PI;
    boltz_unit = 8.61734229648141E-05;  // boltzmann constant ( energy in eV ) 
    time_unit  = 98.2269514139276    ;  // unit of time picosecond => angstrom * ( atomicmassunit / eV ) ** 0.5

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
    time_unit = 1.0;
}
