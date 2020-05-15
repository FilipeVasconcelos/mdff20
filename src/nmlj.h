#ifndef NMLJ_H
#define NMLJ_H
#define NTYPEMAX 16

bool   lsymmetric                   ;
double qlj     [NTYPEMAX][NTYPEMAX] ;
double plj     [NTYPEMAX][NTYPEMAX] ;
double epslj   [NTYPEMAX][NTYPEMAX] ;
double sigmalj [NTYPEMAX][NTYPEMAX] ;
double rcutsq                       ;
double epsp    [NTYPEMAX][NTYPEMAX] ;
double epslj   [NTYPEMAX][NTYPEMAX] ;
double sigsq   [NTYPEMAX][NTYPEMAX] ;
double fc      [NTYPEMAX][NTYPEMAX] ;

double qtwo    [NTYPEMAX][NTYPEMAX] ;
double ptwo    [NTYPEMAX][NTYPEMAX] ;

int  trunctype                         ;  
char trunclabel[3][MAX_LEN+1]       ;

void init_nmlj()                    ;
void info_nmlj()                    ;
void engforce_nmlj_pbc();
#endif
