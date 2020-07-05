#include <math.h>
#include "tt_damp.h"

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
    //printf("ff %e E_TT[order] %e \n",ff,E_TT[order]);
    for (int k=order-1;k>-1;k--){
        ff = ff * br + E_TT[k];
    //    printf("k %d ff %e E_TT %e br %e\n",k,ff,E_TT[k],br);
    }
    //printf("inside TT_damping_functions %e %e %e %e %e %e\n",b,c,r,ff,br,expbdr);
    *f = 1.0 - ff * expbdr;
    //derivative of f checked on sage 18/02/14 worksheet BMHFTD
    *fd = E_TT[order] * pow(br,order) * expbdr * b;
                                                                    
}

void get_TT_damp(){
    E_TT[0] = 1.0;
    for(int i=1;i<(MAX_TT_EXPANSION);i++){
        E_TT[i] = E_TT[i-1] / ( (double) i );
    }
    /*
    for (int i=0;i<MAX_TT_EXPANSION;i++){
        printf("i %e E_TT\n",i,E_TT[i]);
    }
    */
}
