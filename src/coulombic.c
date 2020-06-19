#include <stdlib.h>
#include <stdio.h>
#include "io.h"
#include "field.h"
#include "config.h"
#include "coulombic.h"
#include "kspace.h"
#include "pim.h"

#define DEBUG_COUL
/******************************************************************************/
void allocate_coulombic(){
#ifdef DEBUG_COUL
    printf("inside allocate_coulombic\n");
#endif
    qia=malloc(nion*sizeof(*qia));
    dipia=malloc(nion*sizeof(*dipia));
    quadia=malloc(nion*sizeof(*quadia));
    if (lpim) {
        polia     = malloc(nion*sizeof(*polia));
        invepolia = malloc(nion*sizeof(*invepolia));
        dipia_ind = malloc(nion*sizeof(*dipia_ind));
    }
}
/******************************************************************************/
void free_coulombic(){
    free(qia);
    free(dipia);
    free(quadia);
    free_kspace();
    if (lpim) {
        free(polia);
        free(invepolia);
        free(dipia_ind);
    }
}
/******************************************************************************/
void get_dipoles(double (*mu)[3],double *upol){

    printf("inside get_dipoles %d\n",algo_pim);
    double (*dipia_ind)[3];
    dipia_ind=malloc(nion*sizeof(*dipia_ind));
    switch (algo_pim) {
        case 0 :
            momentpolaSCF(dipia_ind,&upol);
           break; 
        case 1 :
           break; 
    }
    for (int ia=0;ia<nion;ia++){
        for (int i=0;i<3;i++){
            mu[ia][i] = dipia_ind[ia][i];
        }
    }
    free(dipia_ind);
}
/******************************************************************************/
void init_coulombic(){
    allocate_coulombic();
    init_kspace();
    info_coulombic();
}
/******************************************************************************/
void info_coulombic(){
    if (ionode) {
        SEPARATOR;
        printf("coulombic field info\n");
        LSEPARATOR;
        printf("coulombic potential\n");
        LSEPARATOR;
        putchar('\n');
        printf("        qi qj   \n");
        printf(" Vij = -------  \n");
        printf("         rij    \n");
        putchar('\n');
        if (lautoES) {
            printf("ewald summation parameters (from espw) \n");
        }
        else
        {
            printf("ewald summation parameters (from user) \n");
        }
        printf("alpha                            = %f\n",alphaES);
        printf("kmax                             = ");
        for(int k=0;k<3;k++){
            printf("%d ",kES[k]);
        }
        putchar('\n');
        if (!lautoES) {
            printf("relative error (user defined)    : "ee"\n",epsw);
        }
    };
}
/******************************************************************************/
void sample_field_coulombic(double (*ef)[3],double (*efg)[3][3]){

    int mia=10;
    if (nion < mia ) mia=nion;
    if (ionode ) {
        LSEPARATOR;
        printf("              coulombic field sample (first 10)\n");
        LSEPARATOR;
        printf("   ia atype              Ex              Ex              Ez\n");
        for (int ia=0;ia<mia;ia++){
            printf("%5d %5s "ee3"\n",ia,atypia[ia],ef[ia][0],ef[ia][1],ef[ia][2]);
        }
        printf("   ia atype"xx11"EFGxx"xx11"EFGyy"xx11"EFGzz"xx11"EFGxy"xx11"EFGxz"xx11"EFGyz\n");
        for (int ia=0;ia<mia;ia++){
            printf("%5d %5s "ee6"\n",ia,atypia[ia],efg[ia][0][0],efg[ia][1][1],efg[ia][2][2],
                                                   efg[ia][0][1],efg[ia][0][2],efg[ia][0][2]);
        }
        LSEPARATOR;
    }
}

