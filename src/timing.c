#include <stdio.h>
#ifdef MPI
#include <mpi.h>
#endif
#ifdef OMP
#include <omp.h>
#endif
#include <time.h>
#include "field.h"
#include "ewald.h"
#include "timing.h"
#include "io.h"

/******************************************************************************/
void statime(int x){
#if defined(MPI) && !defined(OMP)
    ttt[x]=MPI_Wtime();
#elif defined(OMP)
    ttt[x]=omp_get_wtime();
#else
    ttt[x]=clock();
#endif
}
/******************************************************************************/
void mestime(double *whole,int from,int to){
    *whole+=ttt[from]-ttt[to];
}
/******************************************************************************/
void writime(char* label,int step,int from,int to){
#if defined(MPI) || defined(OMP)
    io_node printf("  step : %-6d  %-6s cpu time : %8.3f (s)\n",step,label,ttt[from]-ttt[to]);
#else
    io_node printf("  step : %-6d  %-6s cpu time : %8.3f (s)\n",step,label,(((double)ttt[from]-ttt[to])/CLOCKS_PER_SEC));
#endif
}
/******************************************************************************/
void writimewhole(char* label, double whole){
#if defined(MPI) || defined(OMP)
    printf("  %-20s  : %12.3f\n",label,whole);
#else
    printf("  %-20s  : %12.3f\n",label,(double) (whole/CLOCKS_PER_SEC));
#endif
}
/******************************************************************************/
void writimewholecount(char* label, double whole, int count){
#if defined(MPI) || defined(OMP)
    printf("  %-20s  : %12.3f  [ %d ]\n",label,whole,count);
#else
    printf("  %-20s  : %12.3f  [ %d ]\n",label,(double) (whole/CLOCKS_PER_SEC),count);
#endif
}

/******************************************************************************/
void info_timing(){

    if (ionode){
        BLANKLINE;
        LSEPARATOR;
        printf("  Timing Info  \n");
        LSEPARATOR;
        writimewhole("Engforce",engforceCPUtime);
        if (lnmlj)   writimewhole("  -> nmlj",engforce_nmljCPUtime);
        if (lbhmft)  writimewhole("  -> bhmft",engforce_bhmftdCPUtime);
        if (lbhmftd) writimewhole("  -> bhmftd",engforce_bhmftdCPUtime);
        if (lcoulombic) {
            if ( lpim ) {
                writimewhole("  -> pim",engforce_getdipolesCPUtime);
                writimewholecount("      -> ewald dir",ewaldDirpimCPUtime,cc_multipole_ES_dir_pim);
                writimewholecount("      -> ewald rec",ewaldRecpimCPUtime,cc_multipole_ES_rec_pim);
            }
            writimewhole("  -> coul",engforce_coulCPUtime);
            writimewholecount("      -> ewald dir",ewaldDirCPUtime,cc_multipole_ES_dir);
            writimewholecount("      -> ewald rec",ewaldRecCPUtime,cc_multipole_ES_rec);
        }
        writimewhole("MD",mdstepCPUtime);
        writimewhole("Verlet List",verletLCPUtime);
        writimewhole("Integration ",propagatorCPUtime);
        writimewhole("Dir<->Kar",dirkardirCPUtime);
        writimewhole("MPI communication",COMMCPUtime);
        writimewhole("test measure",CPUtime);
    }
}
