#include <stdio.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "timing.h"
#include "io.h"

void statime(int x){
#ifdef MPI
    ttt[x]=MPI_Wtime();
#else
    clock_t ttt[x]=clock();
#endif
}
void mestime(double *whole,int from,int to){
    *whole+=ttt[from]-ttt[to];
}
void writime(char* label,int step,int from,int to){
    io_node printf("  step : %d  %s %12.3f\n",step,label,ttt[from]-ttt[to]);
}
void writimewhole(char* label, double whole){
    printf("  %s  : %12.3f\n",label,whole);
}

void info_timing(){

    if (ionode){
        BLANKLINE;
        LSEPARATOR;
        printf("print timing info\n");
        LSEPARATOR;
        writimewhole("MD",mdstepCPUtime);
    }
}
