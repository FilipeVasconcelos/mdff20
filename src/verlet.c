#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "config.h"
#include "verlet.h"
#include "field.h"
#include "io.h"
#include "md.h"
#include "cell.h"
#include "pbc.h"
#include "global.h"
#include "timing.h"
#ifdef DEBUG
    #define DEBUG_VERLET_NB
    #define DEBUG_VERLET_COUL
#endif
//#define DEBUG_
#ifdef DEBUG_
    #define DEBUG_VERLET_NB
    #define DEBUG_VERLET_COUL
#endif
/*******************************************************************************/
VERLETL *allocate_verletlist(char* label){
    if ( ionode ) printf("allocate Verlet List %s\n",label);
    VERLETL *vlist=malloc(sizeof(*vlist));
    vlist->list=malloc(nion*VNLMAX*sizeof(*(vlist->list)));
    vlist->point=malloc((nion+1)*sizeof(*(vlist->point)));
    vlist->label=label;
    for(int k=0;k<nion*VNLMAX;k++){
        vlist->list[k]=0;
    }
    // not the best place to do that
    xvl= malloc(nion*sizeof(*xvl));
    yvl= malloc(nion*sizeof(*yvl));
    zvl= malloc(nion*sizeof(*zvl));
    updatevl = 0;
    return vlist;
}
/*******************************************************************************/
void free_verletlist(char* label){
    free(verlet_nb->list);
    free(verlet_nb->point);
    free(verlet_nb);
    free(verlet_coul->list);
    free(verlet_coul->point);
    free(verlet_coul);
    free(xvl);
    free(yvl);
    free(zvl);
}
/*******************************************************************************/
void gen_pbc_verletlist(){

    double rskinsq_nb=verlet_nb->cut*verlet_nb->cut;
    double rskinsq_coul=verlet_coul->cut*verlet_coul->cut;
    double rxi,ryi,rzi;
    double rxij,ryij,rzij;
    double rijsq;
    int icount_nb,k_nb;
    int icount_coul,k_coul;

    for(int k=0;k<nion*VNLMAX;k++){
        verlet_nb->list[k]=0;
        verlet_coul->list[k]=0;
    }
    for(int k=0;k<nion;k++) {
        verlet_nb->point[k]=0;
        verlet_coul->point[k]=0;
    }

    icount_nb=0;
    icount_coul=0;
    for(int ia=atomDec.iaStart;ia<atomDec.iaEnd;ia++) {
        rxi = rx[ia];
        ryi = ry[ia];
        rzi = rz[ia];
        k_nb=0;
        k_coul=0;
        for(int ja=0;ja<nion;ja++){
            if( ((ia>ja) && ((ia+ja)%2==0)) ||
                ((ia<ja) && ((ia+ja)%2!=0)) ) {
                rxij = rxi - rx[ja];
                ryij = ryi - ry[ja];
                rzij = rzi - rz[ja];
                pbc(&rxij,&ryij,&rzij);
                rijsq = rxij * rxij + ryij * ryij + rzij * rzij;
                // non bonded sphere
                if (rijsq<=rskinsq_nb){
                    verlet_nb->list[icount_nb]=ja;
                    icount_nb+=1;
                    k_nb+=1;
                    if ( (icount_nb < 0) || (icount_nb>=VNLMAX*nion)) {
                        pError("out of bound in gen_verletlist\n");
                    }
                }
                // coulomb sphere
                if (rijsq<=rskinsq_coul){
                    verlet_coul->list[icount_coul]=ja;
                    icount_coul+=1;
                    k_coul+=1;
                    if ( (icount_coul < 0) || (icount_coul>=VNLMAX*nion)) {
                        pError("out of bound in gen_verletlist\n");
                    }
                }
            }
        }
        verlet_nb->point[ia]=icount_nb-k_nb;
        verlet_coul->point[ia]=icount_coul-k_coul;
    }
    verlet_nb->point[atomDec.iaEnd]=icount_nb;
    verlet_coul->point[atomDec.iaEnd]=icount_coul;
#ifdef DEBUG_VERLET_NB
    for(int k=0;k<nion;k++) {printf("verlet list  NB  %d %d\n",k,verlet_nb->list[k]); }
    for(int k=0;k<nion;k++) {printf("verlet point NB  %d %d\n",k,verlet_nb->point[k]); }
#endif
#ifdef DEBUG_VERLET_COUL
    for(int k=0;k<nion;k++) {printf("verlet list  COUL%d %d\n",k,verlet_coul->list[k]); }
    for(int k=0;k<nion;k++) {printf("verlet point COUL%d %d\n",k,verlet_coul->point[k]); }
#endif
}
/*******************************************************************************/
void check_verletlist(){

    double drneimax=0.0,drneimax2=0.0,drnei;
    double rxvl,ryvl,rzvl;
    /***************************************
            cartesian to direct
     ***************************************/
    kardir ( nion , rx , ry , rz , simuCell.B ) ;

    statime(10);
    for(int ia=0;ia<nion;ia++){
        rxvl=rx[ia]-xvl[ia];
        ryvl=ry[ia]-yvl[ia];
        rzvl=rz[ia]-zvl[ia];
        pbc(&rxvl,&ryvl,&rzvl);
        drnei=sqrt(rxvl*rxvl+ryvl*ryvl+rzvl*rzvl);
        if ( drnei > drneimax ) {
            drneimax2=drneimax;
            drneimax=drnei;
        }
        else if ( drnei > drneimax2) {
                drneimax2=drnei;
        }
    }
    if ( drneimax + drneimax2 > skindiff )
    {
        updatevl+=1;
        gen_pbc_verletlist();
        // save coordinates
        for (int ia=0;ia<nion;ia++){
            xvl[ia]=rx[ia];
            yvl[ia]=ry[ia];
            zvl[ia]=rz[ia];
        }
    }
    if (updatevl>1){
        io_pnode printf("  update verlet list : %d/%d T = %f\n",
                           updatevl,istep,((double)istep/updatevl));
    }
    statime(11);
    mestime(&verletLCPUtime,11,10);
    /***************************************
        direct to cartesian
    ***************************************/
    dirkar ( nion , rx , ry , rz , simuCell.A ) ;

}
