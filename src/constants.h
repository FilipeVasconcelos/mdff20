#ifndef CONSTANTS_H
#define CONSTANTS_H
#define NTYPEMAX 16
#define MAX_LEN 80
#define PI (3.14159265358979323846264338327950)

double TPI ;
double boltz_unit;
double time_unit;
double mass   [NTYPEMAX];
int    natmi  [NTYPEMAX];
char*  atypei [NTYPEMAX];  

void gen_constants();
#endif
