#ifndef TIMING_H
#define TIMING_H
#define MAX_INDEX_TIME 32
// 0 -> 1   MD step
// 2 -> 3   engforce
// 4 -> 5   dirkardir
// 6 -> 7   test measure
// 8 -> 9   COMM
// 10 -> 11 verlet
// 12 -> 13 propagator/integration
// 14 -> 15 nmlj
// 16 -> 17 coulomb
// 18 -> 19 ewald dir
// 19 -> 20 ewald rec
// 15 -> 21 bhmftd
#include <time.h>
#if defined(MPI) || defined(OMP)
double ttt[MAX_INDEX_TIME]; /* ttt[0] reference time */
#else
clock_t ttt[MAX_INDEX_TIME]; /* ttt[0] reference time */
#endif
double mdstepCPUtime;
double engforceCPUtime;
double engforce_nmljCPUtime;
double engforce_bhmftdCPUtime;
double engforce_coulCPUtime;
double ewaldDirCPUtime;
double ewaldRecCPUtime;
double dirkardirCPUtime;
double CPUtime; /* tmp measure */
double COMMCPUtime; /* tmp measure */
double verletLCPUtime ;
double propagatorCPUtime;


void info_timing();
void statime(int x);
void mestime(double* whole,int from,int to);
void writime(char* label,int step, int from,int to);
void writimewhole(char* label, double whole);
#endif
