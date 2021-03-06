#ifndef MPI_MDFF_H
#define MPI_MDFF_H
#ifdef MPI
#include <mpi.h>
#endif
int myrank, numprocs;
// Structure atom decomposition parallelism
struct DEC atomDec;

typedef struct DEC
{
    int dimData;
    int iaStart;
    int iaEnd;
    char* label;
}DEC ;

int do_split(int n,int np,int mrank,DEC* dec,char* lab);
void MPI_Allreduce_sumDouble( double *localSum, int ndim);
void MPI_Allreduce_sumDouble_R2(double (*inout)[3], int ndim);
void MPI_Allreduce_sumDouble_R3(double (*inout)[3][3], int ndim);
#endif
