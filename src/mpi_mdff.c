#include <stdio.h>
#include "io.h"
#include "mpi_mdff.h"

int do_split(int n,int np,int mrank,DEC dec,char* lab){
#ifndef MPI
    //quick return
    return 1; 
#endif

    int x,y;
    int istartV[np], iendV[np];
    int splitnumberV[np];

    int imin=1; int imax=n;
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

    dec.iaStart=istartV[mrank];
    dec.iaEnd=iendV[mrank];
    dec.dimData=(iendV[mrank] - istartV[mrank] )+1;
    dec.label=lab;

    io_node printf("paralelisation - %s decomposition\n",dec.label);
    for (int me=0;me<np;me++){
        printf("rank = %d %s %d to %d load : %d \n",me,dec.label,istartV[me],iendV[me],(iendV[me]-istartV[me] +1));
    }



    return 0; 



}
