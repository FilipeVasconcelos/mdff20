#include <stdlib.h>
#include <stdio.h>
#include "io.h"
#include "field.h"
#include "config.h"
#include "multipole.h"
#include "kspace.h"

void allocate_multipole(){
    qia=malloc(nion*sizeof(*qia));
    dipia=malloc(nion*sizeof(*dipia));
    quadia=malloc(nion*sizeof(*quadia));
}

void free_multipole(){
    free(qia);
    free(dipia);
    free(quadia);
}

void init_multipole(){
    allocate_multipole();
    info_multipole();
    init_kspace();
}

void info_multipole(){
    if (ionode) {
        SEPARATOR;
        printf("multipole field info\n");
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
