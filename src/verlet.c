#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "config.h"
#include "verlet.h"
#include "field.h"
#include "io.h"
#include "cell.h"
#include "pbc.h"
#include "global.h"

VERLETL *allocate_verletlist(char* label){
    printf("allocate verletlist\n");
    VERLETL *vlist=malloc(sizeof(*vlist)+ sizeof(nion*VNLMAX*sizeof(int))+ sizeof((nion+1)*sizeof(int)));
    vlist->list=malloc(nion*VNLMAX*sizeof(*(vlist->list)));
    vlist->point=malloc((nion+1)*sizeof(*(vlist->point)));
    vlist->label=label;
    for(int k=0;k<nion*VNLMAX;k++){
        vlist->list[k]=0;
    //    printf("%d %d\n",k,vlist->list[k]);
    }
    // not the best place to do that
    xvl= malloc(nion*sizeof(*xvl));
    yvl= malloc(nion*sizeof(*yvl));
    zvl= malloc(nion*sizeof(*zvl));
    updatevl = 0;
    return vlist;
}
void free_verletlist(char* label){
    free(verlet_nb->list);
    free(verlet_nb->point);
    free(verlet_nb);
    free(xvl);    
    free(yvl);    
    free(zvl);    
}

void gen_pbc_verletlist(){

    double rskinsq_nb=verlet_nb->cut + skindiff;
    rskinsq_nb*=rskinsq_nb;
    double rxi,ryi,rzi;
    double rxij,ryij,rzij;
    double rijsq;
    int icount_nb,k_nb;

    //printf("inside gen_pbc_verletlist %s\n",verlet_nb->label);
    for(int k=0;k<nion*VNLMAX;k++){verlet_nb->list[k]=0; }
    for(int k=0;k<nion;k++) { verlet_nb->point[k]=0; }
    /*************************************** 
            cartesian to direct                    
     ***************************************/
    kardir ( nion , rx , ry , rz , simuCell.B ) ;

    icount_nb=0;
    //printf("%d %d\n",atomDec.iaStart,atomDec.iaEnd);
    for(int ia=atomDec.iaStart;ia<atomDec.iaEnd;ia++) {
        rxi = rx[ia];
        ryi = ry[ia];
        rzi = rz[ia];
        k_nb=0;
        for(int ja=0;ja<nion;ja++){
            if( ((ia>ja) && ((ia+ja)%2==0)) ||
                ((ia<ja) && ((ia+ja)%2!=0)) ) {
                rxij = rxi - rx[ja];
                ryij = ryi - ry[ja];
                rzij = rzi - rz[ja];
                pbc(&rxij,&ryij,&rzij);
                rijsq = rxij * rxij + ryij * ryij + rzij * rzij;
                if (rijsq<=rskinsq_nb){
                //printf("ia %d ja %d %f %f\n",ia,ja,rijsq,rskinsq_nb);
                    verlet_nb->list[icount_nb]=ja;
                    icount_nb+=1;
                    k_nb+=1;
                    if ( (icount_nb < 0) || (icount_nb>=VNLMAX*nion)) {
                        pError("out of bound in gen_verletlist");
                    }
                }
            }
        }
        verlet_nb->point[ia]=icount_nb-k_nb;
    }
    verlet_nb->point[atomDec.iaEnd]=icount_nb;
    /*************************************** 
            direct to cartesian                   
     ***************************************/
    dirkar ( nion , rx , ry , rz , simuCell.A ) ;
    //for(int k=0;k<nion*VNLMAX;k++){printf("verlet list %d %d\n",k,verlet_nb->list[k]); }
    //for(int k=0;k<nion;k++) {printf("verlet point %d %d\n",k,verlet_nb->point[k]); }
}
void check_verletlist(){

    double drneimax=0.0,drneimax2=0.0,drnei;
    double rxvl,ryvl,rzvl;
    /*************************************** 
            cartesian to direct                    
     ***************************************/
    kardir ( nion , rx , ry , rz , simuCell.B ) ;
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
    /*************************************** 
            direct to cartesian                   
     ***************************************/
    dirkar ( nion , rx , ry , rz , simuCell.A ) ;

}
