#ifndef IO_H
#define IO_H
#include <stdbool.h>
#define SEPARATOR for(int i=0;i<61;i++) putchar('=');putchar('\n') 
#define io_node if ( ionode ) 

char *controlfn; 
bool ionode;


void headerstdout(char* time ,int numprocs);
void init_io();
#endif
