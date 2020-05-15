#ifndef CONSTANTS_H
#define CONSTANTS_H
#define NTYPEMAX 16
#define MAX_LEN 80
#define LOWERCASE 97
#define PI        3.1415926535897932385  /* pi  */
#define TPI       6.2831853071795864770  /* 2pi */ 

double radian;

double boltz_unit;
double time_unit;
double press_unit;

double mass   [NTYPEMAX];
int    natmi  [NTYPEMAX];
char*  atypei [NTYPEMAX];  

void gen_constants();
#endif
