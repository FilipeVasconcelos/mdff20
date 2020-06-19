#ifndef FIELD_H
#define FIELD_H
#include <stdbool.h>
#include "constants.h"
bool       lnonbonded; /* non bonded pair potential nmlj bhmftd */
bool       lsymmetric; /* are those non bonded potentials symmetric */
bool            lnmlj; /* n-m lennard-jones potential */
bool           lbhmft; /* Born Huggins-Meyer potential + Fumi-Tossi */
bool          lbhmftd; /* Born Huggins-Meyer potential + Fumi-Tossi + Damping */ 
bool       lcoulombic; /* Coulombic/Electrostatic potential */
bool lpolar[NTYPEMAX]; /* induced moment from polarizability for type */ 
bool             lpim; /* true if one of the above is true */
bool         lpoldamp; /* electric field damping applied to polarizable atoms*/
/* ****************************************************************************/
/*                       EWALD SUMMATION PARAMETERS                           */
/* ****************************************************************************/
bool          lautoES; /* auto-determination of Ewald parameter from epsw ( accuracy) */
double           epsw; /* accuracy of the ewald sum */
double        alphaES; /*Ewald sum */
int            kES[3]; /* kmax of ewald sum in reciprocal space */


/* ****************************************************************************/
/*                              global parameters                             */
/* ****************************************************************************/
double        srcutsq; /* shortrange cut ^2 */
double        lrcutsq; /* longrange cut ^2 */
/* ****************************************************************************/
/*                                prototypes                                  */
/* ****************************************************************************/
bool is_pim();
bool is_damping();
int read_field(char* controlfn);
void init_field(char* controlfn);
void info_field();
void engforce();
#endif
