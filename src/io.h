#ifndef IO_H
#define IO_H
#include <stdbool.h>
#include "color_mdff.h"
#include "constants.h"
#define SEPARATOR  for(int i=0;i<61;i++) putchar('=');putchar('\n') 
#define LSEPARATOR for(int i=0;i<61;i++) putchar('-');putchar('\n') 
#define BSEPARATOR for(int i=0;i<105;i++) putchar('-');putchar('\n') 
#define BLANKLINE  putchar('\n');
#define BSEPF(X) for(int i=0;i<105;i++) fprintf(X,"-");fprintf(X,"\n")
#define io_node if ( ionode ) 
#define io_pnode if ( ionode && iopnode(istep,npas,nprint) ) 
#define pError(X) if (ionode) fprintf(stderr,RED"ERROR : "RES X)

#define OSZHEADER "   step      Time/Temp   Etot/Pressure    Ekin/Pvir_nb  Utot/Pvir_coul      U_nb/Omega     U_coul/Htot\n"
int oszcall;
char *controlfn; 
bool ionode;

void headerstdout(char* time ,int numprocs);
void init_io();
int iopnode(int step,int npas, int nprint);
int check_FTstring(char* label,char buffer[MAX_LEN+1] );
#endif
