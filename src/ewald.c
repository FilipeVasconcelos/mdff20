#include <stdlib.h>
#include "config.h"
#include "tensors.h"
#include "ewald.h"

void get_monopoles(){
}

void get_dipoles(){
}

void get_quadrupoles(){
}


void alloc_multipole(){
}


void multipole_ewald_sum(){

    double u_dir, u_rec;
    double (*ef_dir)[3], (*ef_rec)[3];
    double (*efg_dir)[3][3], (*efg_rec)[3][3];
    double *fx_dir, *fy_dir, *fz_dir;
    double *fx_rec, *fy_rec, *fz_rec;
    double tau_dir[3][3],tau_rec[3][3];

    ef_dir=malloc(nion*sizeof(*ef_dir));
    efg_dir=malloc(nion*sizeof(*efg_dir));
    fx_dir=malloc(nion*sizeof(*fx_dir));
    fy_dir=malloc(nion*sizeof(*fy_dir));
    fz_dir=malloc(nion*sizeof(*fz_dir));
    multipole_ewald_sum_direct(u_dir, ef_dir, efg_dir, fx_dir, fy_dir, fz_dir, tau_dir);

    ef_rec=malloc(nion*sizeof(*ef_rec));
    efg_rec=malloc(nion*sizeof(*efg_rec));
    fx_rec=malloc(nion*sizeof(*fx_rec));
    fy_rec=malloc(nion*sizeof(*fy_rec));
    fz_rec=malloc(nion*sizeof(*fz_rec));
    multipole_ewald_sum_reciprocal(u_rec, ef_rec, efg_rec, fx_rec, fy_rec, fz_rec, tau_rec);

    free(ef_dir);
    free(ef_rec);
    free(efg_dir);
    free(efg_rec);
    free(fx_dir);free(fy_dir);free(fz_dir);
    free(fx_rec);free(fy_rec);free(fz_rec);
}


void multipole_ewald_sum_direct(double u_dir  , double (*ef_dir)[3], double (*efg_dir)[3][3], 
                                double *fx_dir, double *fy_dir, double *fz_dir , double tau_dir[3][3]){

    TENSOR_RK1 T0;
    TENSOR_RK1 T1;
    TENSOR_RK2 T2;
    TENSOR_RK3 T3;
    TENSOR_RK4 T4;

    /*************************************** 
            cartesian to direct                    
     ***************************************/
    kardir ( nion , rx , ry , rz , simuCell.B ) ;

    for(ia=atomDec.iaStart;ia<atomDec.iaEnd;ia++) {
        rxi = rx[ia];
        ryi = ry[ia];
        rzi = rz[ia];
        if (lverletL) {
            jb=verlet_coul->point[ia];
            je=verlet_coul->point[ia+1];
        }
        else{
            jb = 0;
            je = nion;
        }
    }



    /*************************************** 
            direct to cartesian                   
     ***************************************/
    dirkar ( nion , rx , ry , rz , simuCell.A ) ;
}

void multipole_ewald_sum_reciprocal(double u_rec  , double (*ef_rec)[3], double (*efg_rec)[3][3],                                                                 double *fx_rec, double *fy_rec, double *fz_rec , double tau_rec[3][3]){
}

