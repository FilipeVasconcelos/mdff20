#ifndef NMLJ_H
#define NMLJ_H
#define NTYPEMAX 16

double qlj    [NTYPEMAX][NTYPEMAX] ;
double plj    [NTYPEMAX][NTYPEMAX] ;
double epslj  [NTYPEMAX][NTYPEMAX] ;
double sigmalj [NTYPEMAX][NTYPEMAX] ;
//double rcutsq [NTYPEMAX][NTYPEMAX] ;
double rcutsq ;
double epsp   [NTYPEMAX][NTYPEMAX] ;
double epslj   [NTYPEMAX][NTYPEMAX] ;
double sigsq  [NTYPEMAX][NTYPEMAX] ;
double fc     [NTYPEMAX][NTYPEMAX] ;

double qtwo    [NTYPEMAX][NTYPEMAX] ;
double ptwo   [NTYPEMAX][NTYPEMAX] ;

void init_nmlj();
void engforce_nmlj_pbc();
#endif
