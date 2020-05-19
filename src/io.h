#ifndef IO_H
#define IO_H
#include <stdbool.h>
#include "color_mdff.h"
#include "constants.h"
#define SEPARATOR  for(int i=0;i<61;i++) putchar('=');putchar('\n') 
#define LSEPARATOR for(int i=0;i<61;i++) putchar('-');putchar('\n') 
#define BLANKLINE  putchar('\n');
#define io_node if ( ionode ) 
#define io_pnode if ( (ionode) && (istep % nprint == 0) ) 
#define pError(X) if (ionode) fprintf(stderr,RED"ERROR : "RES X)

char *controlfn; 
bool ionode;

void headerstdout(char* time ,int numprocs);
void init_io();
int check_FTstring(char* label,char buffer[MAX_LEN+1] );
#endif
