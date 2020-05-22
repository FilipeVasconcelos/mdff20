#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "constants.h"

/******************************************************************************/
void gen_constants()
{

    boltz_unit = 8.61734229648141E-05F;  // boltzmann constant ( energy in eV ) 
    time_unit  = 98.2269514139276F    ;  // unit of time picosecond => angstrom * ( atomicmassunit / eV ) ** 0.5
    press_unit = 0.00624150964712042F ;  // GPa => internal unit of pressure (eV/A^3)
    rho_unit   = 1.660538782F         ;

    radian   = 180. / PI;     
    onethird = 1./3.0;

    strcpy(allwd_FT_str[0 ],"yes");
    strcpy(allwd_FT_str[1 ],"no");
    strcpy(allwd_FT_str[2 ],"y");
    strcpy(allwd_FT_str[3 ],"n");
    strcpy(allwd_FT_str[4 ],"Y");
    strcpy(allwd_FT_str[5 ],"N");
    strcpy(allwd_FT_str[6 ],"Yes");
    strcpy(allwd_FT_str[7 ],"No");
    strcpy(allwd_FT_str[8 ],"true");
    strcpy(allwd_FT_str[9 ],"false");
    strcpy(allwd_FT_str[10],"t");
    strcpy(allwd_FT_str[11],"f");
    strcpy(allwd_FT_str[12],"T");
    strcpy(allwd_FT_str[13],"F");
    strcpy(allwd_FT_str[14],"True");
    strcpy(allwd_FT_str[15],"False");
}
/******************************************************************************/
void reduced_units()
{
    boltz_unit = 1.0;
    time_unit  = 1.0;
    press_unit = 1.0;
    rho_unit   = 1.0;
}
