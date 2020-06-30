#ifndef TTDAMP_H
#define TTDAMP_H

#define MAX_TT_EXPANSION 9
double E_TT[MAX_TT_EXPANSION];

void TT_damping_functions(double b,double c,double r,double *f,double *fd,int order);
void get_TT_damp();
#endif
