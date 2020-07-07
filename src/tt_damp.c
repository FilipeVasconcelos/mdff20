#include <math.h>
#include "tt_damp.h"
#ifdef DEBUG
    #define DEGUG_E_TT
#endif
//#define DEBUG_
#ifdef DEBUG_
    #define DEGUG_E_TT
#endif
/******************************************************************************/
/* Tang-Toennies damping function                                             */
/* fmv adapted from fortran code (Feb 2014)                                   */
/* f : damping function                                                       */
/* fd: first derivative of f                                                  */
/******************************************************************************/
void TT_damping_functions(double b,double c,double r,double *f,double *fd,int order){

    double expbdr, br;
    double ff;

    br = b * r;
    expbdr = exp(-br) * c;

    ff = E_TT[order];
    for (int k=order-1;k>-1;k--){
        ff = ff * br + E_TT[k];
    }
    *f = 1.0 - ff * expbdr;
    //derivative of f checked on sage 18/02/14 worksheet BMHFTD
    *fd = E_TT[order] * pow(br,order) * expbdr * b;

}

void get_TT_damp(){
    E_TT[0] = 1.0;
    for(int i=1;i<(MAX_TT_EXPANSION);i++){
        E_TT[i] = E_TT[i-1] / ( (double) i );
    }
#ifdef DEBUG_E_TT
    for (int i=0;i<MAX_TT_EXPANSION;i++){
        printf("i %e E_TT\n",i,E_TT[i]);
    }
#endif
}
