#ifndef COULOMBIC_H
#define COULOMBIC_H
#include <stdbool.h>
bool lqch;
bool ldip;
bool lqua;

/* ****************************************************************************/
/*                                prototypes                                  */
/* ****************************************************************************/
void sample_field_coulombic(double (*ef)[3],double (*efg)[3][3]);
void init_coulombic();
void info_coulombic();
void free_coulombic();
void get_dipoles(double (*mu)[3],double *upol);
#endif
