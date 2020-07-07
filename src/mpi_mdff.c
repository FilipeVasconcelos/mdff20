#include <stdio.h>
#include <stdlib.h>
#include "io.h"
#include "mpi_mdff.h"

/******************************************************************************/
int do_split(int n,int np,int mrank, DEC* dec, char* lab){

    int x,y;
    int istartV[np], iendV[np];
    int splitnumberV[np];
    int imin=0; int imax=n-1;
    int isteps= (imax-imin)+1;
    x=isteps/np;
    y=isteps%np;

    for (int me=0;me<np;me++) {
        if ((me==0) || (me>y)) {
            splitnumberV[me]=x;
        }
        else if ( (me>0) || (me<(y+1))) {
            splitnumberV[me]=x+1;
        }
    }
    for (int me=0;me<np;me++){
        if (me == 0) {
            istartV[me]=imin;
            iendV[me]=imin+x-1;
        }
        else if (me>0) {
            istartV[me]=iendV[me-1]+1;
            iendV[me]=istartV[me]+splitnumberV[me]-1;
        }
    }
    dec->iaStart=istartV[mrank];
    dec->iaEnd=iendV[mrank]+1;
    dec->dimData=(iendV[mrank] - istartV[mrank] )+1;
    dec->label=lab;
    if (ionode) {
        SEPARATOR;
        printf("paralelisation info\n");
        LSEPARATOR;
        printf("%s decomposition\n",dec->label);
        for (int me=0;me<np;me++){
            printf("rank = %d %s from %4d to %4d load : %d \n",\
            me,dec->label,istartV[me],iendV[me],(iendV[me]-istartV[me] +1));
        }
        putchar('\n');
    }
    return 0;
}

/******************************************************************************/
void MPI_Allreduce_sumDouble( double *localSum, int ndim){
    if (ndim ==1 ) {
        double globalSum;
#ifdef MPI
        MPI_Allreduce(localSum, &globalSum, ndim, MPI_DOUBLE, MPI_SUM ,MPI_COMM_WORLD);
#else
        globalSum=*localSum;
#endif
        *localSum=globalSum;
    }
    else
    {
        double *globalSum=malloc(ndim*sizeof(*globalSum));
#ifdef MPI
        MPI_Allreduce(localSum, globalSum, ndim, MPI_DOUBLE, MPI_SUM ,MPI_COMM_WORLD);
#else
        for (int i=0;i<ndim;i++){globalSum[i]=localSum[i];}
#endif
        for (int i=0;i<ndim;i++){localSum[i]=globalSum[i];}
        free(globalSum);
    }
}


//void MPI_Allreduce_sumArrayDouble( double *localSum, double *globalSum, int ndim ) {
//    MPI_Allreduce(localsum,globalSum,)
//}



