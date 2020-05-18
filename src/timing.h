#ifndef TIMING_H
#define TIMING_H

#define MAX_INDEX_TIME 4 

#ifdef MPI
double ttt[MAX_INDEX_TIME]; /* ttt[0] reference time */
#else
clock_t ttt[MAX_INDEX_TIME]; /* ttt[0] reference time */
#endif
double mdstepCPUtime;


void info_timing();
void statime(int x);
void mestime(double* whole,int from,int to);
void writime(char* label,int step, int from,int to);
void writimewhole(char* label, double whole);
#endif
