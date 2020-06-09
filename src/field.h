#ifndef FIELD_H
#define FIELD_H
#include <stdbool.h>
bool lnonbonded; /* non bonded pair potential nmlj bhmftd */
bool lsymmetric; /* are those non bonded potentials symmetric */
bool      lnmlj; /* n-m lennard-jones potential */
bool     lbhmft; /* Born Huggins-Meyer potential + Fumi-Tossi */
bool    lbhmftd; /* Born Huggins-Meyer potential + Fumi-Tossi + Damping */ 
bool lcoulombic; /* Coulombic/Electrostatic potential */


bool    lautoES; /* auto-determination of Ewald parameter from epsw ( accuracy) */
double     epsw; /* accuracy of the ewald sum */
double  alphaES; /*Ewald sum */
int      kES[3]; /* kmax of ewald sum in reciprocal space */


int read_field(char* controlfn);
void init_field(char* controlfn);
void info_field();
void engforce();
#endif
