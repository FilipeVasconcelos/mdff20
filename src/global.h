#ifndef GLOBAL_H
#define GLOBAL_H
int           Fposff;
#include <stdbool.h>
bool        lverletL;
bool         lstatic;
double cutshortrange;
double  cutlongrange;
double      skindiff;
/* ****************************************************************************/
/*                                prototypes                                  */
/* ****************************************************************************/
void init_global(char* controlfn);
#endif
