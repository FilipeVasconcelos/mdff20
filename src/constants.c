#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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

    strcpy(allowed_Fstring[0],"no");
    strcpy(allowed_Fstring[1],"No");
    strcpy(allowed_Fstring[2],"n");
    strcpy(allowed_Fstring[3],"N");
    strcpy(allowed_Fstring[4],"F");
    strcpy(allowed_Fstring[5],"f");
    strcpy(allowed_Fstring[6],"False");
    strcpy(allowed_Fstring[7],"false");
    strcpy(allowed_Tstring[0],"yes");
    strcpy(allowed_Tstring[1],"Yes");
    strcpy(allowed_Tstring[2],"y");
    strcpy(allowed_Tstring[3],"Y");
    strcpy(allowed_Tstring[4],"T");
    strcpy(allowed_Tstring[5],"t");
    strcpy(allowed_Tstring[6],"True");
    strcpy(allowed_Tstring[7],"true");
}


void reduced_units()
{
    boltz_unit = 1.0;
    time_unit  = 1.0;
    press_unit = 1.0;
    rho_unit   = 1.0;
}
