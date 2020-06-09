#ifndef BHMFTD_H
#define BHMFTDH
#define NTYPEMAX 16
#include <stdbool.h>
/* -------------------------------*/
// from input 
double Abhmftd  [NTYPEMAX][NTYPEMAX] ;
double Bbhmftd  [NTYPEMAX][NTYPEMAX] ;
double Cbhmftd  [NTYPEMAX][NTYPEMAX] ;
double Dbhmftd  [NTYPEMAX][NTYPEMAX] ;
double BDbhmftd [NTYPEMAX][NTYPEMAX] ;
/* -------------------------------*/
// from init
// global 
int    trunctyp                      ;  
double utail                         ; /* long range correction to energy   */
double ptail                         ; /* long range correction to pressure */

/* function prototypes */
void init_bhmftd(char* controlfn)    ;
void info_bhmftd()                   ;
void engforce_bhmftd_pbc(double *u, double *pvir, double tau[3][3])            ;
#endif
