#ifndef FIELD_H
#define FIELD_H
#include <stdbool.h>
bool lnonbonded;
bool      lnmlj;
bool      lcoul;


bool    lautoES; /* auto-determination of Ewald parameter from epsw ( accuracy) */
double     epsw; /* accuracy of the ewald sum */
double  alphaES; /*Ewald sum */
int   kES[3]; /* kmax of ewald sum in reciprocal space */


int read_field(char* controlfn);
void init_field(char* controlfn);
void info_field();
void engforce();
#endif
