#ifndef NMLJ_H
#define NMLJ_H
#define NTYPEMAX 16
#include <stdbool.h>
bool   lsymmetric                   ;
/* -------------------------------*/
// from input 
double qlj     [NTYPEMAX][NTYPEMAX] ;
double plj     [NTYPEMAX][NTYPEMAX] ;
double epslj   [NTYPEMAX][NTYPEMAX] ;
double sigmalj [NTYPEMAX][NTYPEMAX] ;
char trunclabel[3][MAX_LEN+1]       ;
/* -------------------------------*/
// from init
double sigsq   [NTYPEMAX][NTYPEMAX] ;
double epsqp   [NTYPEMAX][NTYPEMAX] ;
double fc      [NTYPEMAX][NTYPEMAX] ;
double qtwo    [NTYPEMAX][NTYPEMAX] ;
double ptwo    [NTYPEMAX][NTYPEMAX] ;
double uc      [NTYPEMAX][NTYPEMAX] ;
double c1      [NTYPEMAX][NTYPEMAX] ;
double c2      [NTYPEMAX][NTYPEMAX] ;
// global 
double rcutsq                       ; /* cutshortrange^2 */
int  trunctype                      ;  
double utail                        ; /* long range correction to energy   */
double ptail                        ; /* long range correction to pressure */

/* function prototypes */
void init_nmlj()                    ;
void info_nmlj()                    ;
void engforce_nmlj_pbc(double *u, double *pvir, double tau[3][3])            ;
#endif
