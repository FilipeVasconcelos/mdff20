#ifndef FIELD_H
#define FIELD_H
#include <stdbool.h>
bool lnmlj;
double cutshortrange;
double skindiff;



int read_field(char* controlfn);
void init_field(char* controlfn);
void info_field();
void engforce();
#endif
