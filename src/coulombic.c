#include <stdlib.h>
#include <stdio.h>
#include "io.h"
#include "field.h"
#include "config.h"
#include "coulombic.h"
#include "kspace.h"
#include "pim.h"

#ifdef DEBUG
#define DEBUG_COUL
#define DEBUG_DIPOLE
#endif

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
    for(int ia=0;ia<nion;ia++){
        qia[ia]=0.0;
        for(int i=0;i<3;i++){
            dipia[ia][i]=0.0;
            for(int j=0;j<3;j++){
                quadia[ia][i][j]=0.0;
            }
            if (lpim){
                for(int i=0;i<3;i++){
                    for(int k=0;k<MAX_EXTRAPOLATE_ORDER+1;k++){
                        dipia_ind[ia][k][i]=0.0;
                    }
                    for(int j=0;j<3;j++){
                        polia[ia][i][j]=0.0;
                        invepolia[ia][i][j]=0.0;
                    }
                }
            }
        }
    }
}
/******************************************************************************/
void free_coulombic(){
    free(qia);
    free(dipia);
    free(quadia);
    free_kspace();
    if (lpim) {
        printf("LPIM \n");
        free(polia);
        free(invepolia);
        free(dipia_ind);
    }
}
/******************************************************************************/
void get_dipoles(double (*mu)[3],double *upol){

#ifdef DEBUG_DIPOLE
    printf("inside get_dipoles\n");
#endif
    // save current forces
    double *fx_s,*fy_s,*fz_s;
    fx_s=malloc(nion*sizeof(*fx_s));
    fy_s=malloc(nion*sizeof(*fy_s));
    fz_s=malloc(nion*sizeof(*fz_s));
    for (int ia=0;ia<nion;ia++){
        fx_s[ia]=fx[ia];
        fy_s[ia]=fy[ia];
        fz_s[ia]=fz[ia];
    }
    // static dipoles
    for (int ia=0;ia<nion;ia++){
        for (int i=0;i<3;i++){
            mu[ia][i] = dipia[ia][i];
        }
    }
    // induced dipoles from polarizabilities
    if (lpol) {
        // allocate and set to zero
        double (*mu_ind)[3];
        mu_ind=malloc(nion*sizeof(*mu_ind));
        for (int ia=0;ia<nion;ia++){
            for (int i=0;i<3;i++){
               mu_ind[ia][i]=0.0;
            }
        }
        //pim algo !
        switch (algo_pim) {
            case 0 :
                momentpolaSCF(mu_ind,upol);
               break; 
            case 1 :
               break; 
        }
        for (int ia=0;ia<nion;ia++){
            for (int i=0;i<3;i++){
                mu[ia][i] += mu_ind[ia][i]; 
            }
#ifdef DEBUG_DIPOLE
            printf("ia %d mu_ind %e %e %e\n",ia,mu_ind[ia][0],mu_ind[ia][1],mu_ind[ia][2]);
#endif
        }
        free(mu_ind);
    }
    //restore current forces
    for (int ia=0;ia<nion;ia++){
        fx[ia]=fx_s[ia];
        fy[ia]=fy_s[ia];
        fz[ia]=fz_s[ia];
    }
    free(fx_s);
    free(fy_s);
    free(fz_s);
#ifdef DEBUG_DIPOLE
    printf("out of get_dipoles\n");
#endif
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
        if ( (ldip) || (lqch) ) {
            if (lqch) {
                printf("charge-charge interaction : \n");
                putchar('\n');
                printf("        qi qj   \n");
                printf(" Vij = -------  \n");
                printf("         rij    \n");
                putchar('\n');
            }
            if (ldip) {
                printf("dipole-dipole interaction : \n");
                putchar('\n');
                printf("                /                  -->  -->     -->   -->   \\\n"); 
                printf("           1    | -->   -->      ( mui .rij ) ( rij . muj )  |\n");
                printf(" Vij =  ------- | mui . muj - 3 ---------------------------  |\n");
                printf("         rij^3  \\                            rij^2          /\n"); 
                putchar('\n');
            }
            if ( (ldip) && (lqch) ) {
                printf("charge-dipole interaction : \n");
                putchar('\n');
                printf("             -->   -->   \n");  
                printf("          qi muj . rij   \n"); 
                printf(" Vij =  ---------------- \n"); 
                printf("             rij^3       \n");  
                putchar('\n');
            }
        }
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
        putchar('\n');
    }
}

