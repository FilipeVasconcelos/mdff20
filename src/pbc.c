#include "cell.h"
#include <math.h>
#include "math_mdff.h"

void pbc(double *rxij,double *ryij,double *rzij){
    double sxij = *rxij - nint ( *rxij );
    double syij = *ryij - nint ( *ryij );
    double szij = *rzij - nint ( *rzij );
    *rxij = sxij * simuCell.A[0][0] + syij * simuCell.A[1][0] + szij * simuCell.A[2][0];
    *ryij = sxij * simuCell.A[0][1] + syij * simuCell.A[1][1] + szij * simuCell.A[2][1];
    *rzij = sxij * simuCell.A[0][2] + syij * simuCell.A[1][2] + szij * simuCell.A[2][2];
}
