#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>
#include "io.h"

void invertmatrix3x3(double A[3][3]){
    int n = 3;
    int nn= n*n;
//    int pivotarr[3];
    int ierr;
    double WORK[9];
    int ipiv[3];

    dgetrf_(&n,&n,A[0],&n,ipiv,&ierr);
    if (ierr<0) {
        pError("call to DGETRF failed in invertmatrix3x3\n");
    }
    dgetri_(&n,A[0],&n,ipiv,WORK,&nn,&ierr);
    if (ierr<0) {
        pError("call to DGETRI failed in invertmatrix3x3\n");
    }

}
