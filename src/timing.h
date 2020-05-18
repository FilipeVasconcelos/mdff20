#ifndef TIMING_H
#define TIMING_H

#define MAX_INDEX_TIME 16 
// 0 -> 1 MD step
// 2 -> 3 engforce
// 4 -> 5 dirkardir
// 6 -> 7 test measure
// 8 -> 9 COMM

#ifdef MPI
double ttt[MAX_INDEX_TIME]; /* ttt[0] reference time */
#else
clock_t ttt[MAX_INDEX_TIME]; /* ttt[0] reference time */
#endif
double mdstepCPUtime;
double engforceCPUtime;
double dirkardirCPUtime;
double CPUtime; /* tmp measure */
double COMMCPUtime; /* tmp measure */


void info_timing();
void statime(int x);
void mestime(double* whole,int from,int to);
void writime(char* label,int step, int from,int to);
void writimewhole(char* label, double whole);
#endif
