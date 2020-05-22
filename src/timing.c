#include <stdio.h>
#ifdef MPI
#include <mpi.h>
#endif
#ifdef OMP 
#include <omp.h>
#endif
#include <time.h>
#include "timing.h"
#include "io.h"

void statime(int x){
#ifdef MPI
    ttt[x]=MPI_Wtime();
#elif OMP
//    printf("omp_get_wtime\n");
    ttt[x]=omp_get_wtime();
#else
    ttt[x]=clock();
#endif
}
void mestime(double *whole,int from,int to){
    *whole+=ttt[from]-ttt[to];
}
void writime(char* label,int step,int from,int to){
#if defined(MPI) || defined(OMP)
    io_node printf("  step : %-6d  %-6s cpu time : %8.3f (s)\n",step,label,ttt[from]-ttt[to]);
#else
    io_node printf("  step : %-6d  %-6s cpu time : %8.3f (s)\n",step,label,(((double)ttt[from]-ttt[to])/CLOCKS_PER_SEC));
#endif
}
void writimewhole(char* label, double whole){
#if defined(MPI) || defined(OMP)
    printf("  %-20s  : %12.3f\n",label,whole);
#else
    printf("  %-20s  : %12.3f\n",label,(double) (whole/CLOCKS_PER_SEC));
#endif
}

void info_timing(){

    if (ionode){
        BLANKLINE;
        LSEPARATOR;
        printf("  Timing Info  \n");
        LSEPARATOR;
        writimewhole("MD",mdstepCPUtime);
        writimewhole("Engforce",engforceCPUtime);
        writimewhole("Verlet List",verletLCPUtime);
        writimewhole("Integration ",propagatorCPUtime);
        writimewhole("Dir<->Kar",dirkardirCPUtime);
        writimewhole("MPI communication",COMMCPUtime);
        writimewhole("test measure",CPUtime);
    }
}
