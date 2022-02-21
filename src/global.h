#ifndef GLOBAL_H
#define GLOBAL_H
#include <stdbool.h>
#include "constants.h"
int           Fposff;  /* Format of POSFF                                        [default : 0 ]    */
bool        lverletL;  /* use verlet list                                        [default : true]  */
double      skindiff;  /* skin distance in verlet list */
bool         lstatic;  /* perform static  calculation                            [default : false] */
bool        lreduced;  /* print out quantities in reduced units see constant.h */
bool            lrdf;  /* calcul radial distribution function */
double cutshortrange;  /* cut off for shortrange interaction (nmlj,bhmftd) */
double  cutlongrange;  /* cutl off for longrange interaction (coulombic) */
bool        lpstress;  /* print stress tensor for short and longrange inteaction [default : false] */

#define ALLWD_FORMAT_POSFF_STR 3
char allwd_Fposff[ALLWD_FORMAT_POSFF_STR][MAX_LEN+1];


/* ****************************************************************************/
/*                                prototypes                                  */
/* ****************************************************************************/
void init_global(char* controlfn);
#endif
