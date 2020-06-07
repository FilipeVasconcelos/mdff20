#include <stdlib.h>
#include <stdio.h>
#include "io.h"
#include "field.h"
#include "config.h"
#include "coulombic.h"
#include "kspace.h"

void sample_field_coulombic(double (*ef)[3],double (*efg)[3][3]){

    int mia=10;
    if (nion < mia ) mia=nion;
    if (ionode ) {
        LSEPARATOR;
        printf("              coulombic fiedl sample (first 10)\n");
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


void allocate_coulombic(){
    qia=malloc(nion*sizeof(*qia));
    dipia=malloc(nion*sizeof(*dipia));
    quadia=malloc(nion*sizeof(*quadia));
}

void free_coulombic(){
    free(qia);
    free(dipia);
    free(quadia);
}

void init_coulombic(){
    allocate_coulombic();
    info_coulombic();
    init_kspace();
}

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
        printf("ewald summation parameters\n");
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
