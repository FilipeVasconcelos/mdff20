#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "field.h"
#include "kspace.h"
#include "cell.h"
#include "tools.h"

/******************************************************************************/
void init_kspace(){

    int nk; /* number k points */
    strcpy(kcoul.meshlabel,"kpt-coul");
    for(int k=0;k<3;k++){
        kcoul.kmax[k]=kES[k];
    }
    nk = kcoul.kmax[2] + kcoul.kmax[1] * (2*kcoul.kmax[2]+1) + kcoul.kmax[0] * (2*kcoul.kmax[1]+1) * (2*kcoul.kmax[2]+1);
    kcoul.nk = nk;
    printf("number of k-point %d\n",nk);
    kcoul.kx=malloc(nk*sizeof(*kcoul.kx));
    kcoul.ky=malloc(nk*sizeof(*kcoul.ky));
    kcoul.kz=malloc(nk*sizeof(*kcoul.kz));
    kcoul.kk=malloc(nk*sizeof(*kcoul.kk));
    kcoul.Ak=malloc(nk*sizeof(*kcoul.Ak));
    kcoul.kcoe=malloc(nk*sizeof(*kcoul.kcoe));
    set_param_kmesh(&kcoul,alphaES);
    reorder_kmesh(&kcoul);
}

/******************************************************************************/
void free_kspace(){
    free(kcoul.kx);
    free(kcoul.ky);
    free(kcoul.kz);
    free(kcoul.kk);
    free(kcoul.Ak);
    free(kcoul.kcoe);
}

/******************************************************************************/
void set_param_kmesh(KMESH *km,double alpha){

    int nk;
    int nymin,nzmin;
    double kk,kx,ky,kz;
    double alpha2=alpha*alpha;
    nk=0;
    for (int nx=0;nx<km->kmax[0]+1;nx++){
        if (nx == 0) {
            nymin = 0;
        }
        else {
            nymin = - km->kmax[1];
        }
        for (int ny = nymin;ny<km->kmax[1]+1;ny++){
            if ( (nx == 0) && (ny == 0)) {
                nzmin = 1;
            }
            else {
                nzmin = - km->kmax[2];
            }
            for (int nz = nzmin;nz<km->kmax[2]+1;nz++){

                kx = TPI * ((double) nx*simuCell.B[0][0] + (double) ny*simuCell.B[0][1] + (double) nz*simuCell.B[0][2]);
                ky = TPI * ((double) nx*simuCell.B[1][0] + (double) ny*simuCell.B[1][1] + (double) nz*simuCell.B[1][2]);
                kz = TPI * ((double) nx*simuCell.B[2][0] + (double) ny*simuCell.B[2][1] + (double) nz*simuCell.B[2][2]);
                kk = kx*kx + ky*ky + kz*kz;
                km->kx[nk]=kx;
                km->ky[nk]=ky;
                km->kz[nk]=kz;
                km->kk[nk]=kk;
                km->Ak[nk]= exp(-kk*0.25/alpha2)/kk;
                km->kcoe[nk]= 2.0 * ( 1.0 / kk + 1.0 / alpha2 / 4.0 );
                nk = nk + 1;
            }
        }
    }
    km->nk=nk;
}


/******************************************************************************/
void reorder_kmesh(KMESH *km){

    int *index;
    int is;
    int nk=km->nk; /* number of k-points */
    double *tAk, *tkcoe, *tkx, *tky, *tkz, *tkk;
    tkx=malloc(nk*sizeof(*tkx));
    tky=malloc(nk*sizeof(*tky));
    tkz=malloc(nk*sizeof(*tkz));
    tkk=malloc(nk*sizeof(*tkk));
    tAk=malloc(nk*sizeof(*tAk));
    tkcoe=malloc(nk*sizeof(*tkcoe));

    /* sort index array according to k**2 module */
    index=malloc(nk*sizeof(*index));
    mainsort(km->kk, nk, index);

    /* save before sort */
    for (int ik=0;ik<nk;ik++){
        tkx[ik]=km->kx[ik];
        tky[ik]=km->ky[ik];
        tkz[ik]=km->kz[ik];
        tkk[ik]=km->kk[ik];
        tAk[ik]=km->Ak[ik];
        tkcoe[ik]=km->kcoe[ik];
    }
    /* using sorted index */
    for (int ik=0;ik<nk;ik++){
        is=index[ik];
        km->kx[ik]=tkx[is];
        km->ky[ik]=tky[is];
        km->kz[ik]=tkz[is];
        km->kk[ik]=tkk[is];
        km->Ak[ik]=tAk[is];
        km->kcoe[ik]=tkcoe[is];
    }

    free(tkx);
    free(tky);
    free(tkz);
    free(tkk);
    free(tAk);
    free(tkcoe);

}
