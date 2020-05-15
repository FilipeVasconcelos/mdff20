#ifndef IO_H
#define IO_H
#include <stdbool.h>


#define SEPARATOR  for(int i=0;i<61;i++) putchar('=');putchar('\n') 
#define LSEPARATOR for(int i=0;i<61;i++) putchar('-');putchar('\n') 
#define BLANKLINE  putchar('\n');

#define io_node if ( ionode ) 

#define pError(X) if (ionode) fprintf(stderr,RED"ERROR : "RES X)


char *controlfn; 
bool ionode;


void headerstdout(char* time ,int numprocs);
void init_io();
#endif
