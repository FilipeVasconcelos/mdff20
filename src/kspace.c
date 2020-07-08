#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "field.h"
#include "kspace.h"
#include "cell.h"
#include "tools.h"
#include "config.h"

//#define DEBUG_STRUCT_FACT

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
    kcoul.rhon_R=malloc(nk*sizeof(*kcoul.rhon_R));
    kcoul.rhon_I=malloc(nk*sizeof(*kcoul.rhon_I));
    kcoul.ckria=malloc(nk*sizeof(**kcoul.ckria));
    kcoul.skria=malloc(nk*sizeof(**kcoul.skria));
    kcoul.str=malloc(nk*sizeof(*kcoul.str));
    for (int ik=0;ik< kcoul.nk;ik++) {
        kcoul.ckria[ik]=malloc(nion*sizeof(*kcoul.ckria));
        kcoul.skria[ik]=malloc(nion*sizeof(*kcoul.skria));
    }
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
    free(kcoul.rhon_R);
    free(kcoul.rhon_I);
    free(kcoul.str);
    for (int ik=0;ik<kcoul.nk;ik++) {
        free(kcoul.ckria[ik]);
        free(kcoul.skria[ik]);
    }
    free(kcoul.ckria);
    free(kcoul.skria);

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


/******************************************************************************/
/* Multipole expansion of the Ewald sum  in reciprocal space                  */
/******************************************************************************/
void struct_fact(){

#ifdef DEBUG_STRUCT_FACT
    printf("inside struct_fact\n");
#endif
    double kx,ky,kz,Ak,kcoe;
    double rxi,ryi,rzi;
    double k_dot_r;

    int ik,ia;

    #pragma omp parallel default(none) \
                         shared(rx,ry,rz,kcoul,nion) \
                         private(ik,ia,kx,ky,kz,Ak,kcoe,rxi,ryi,rzi,k_dot_r)
    {
        for (ik=kcoul.kptDec.iaStart;ik<kcoul.kptDec.iaEnd;ik++){

#ifdef DEBUG_STRUCT_FACT_
            printf("k %d %d %d\n",ik,kcoul.kptDec.iaStart,kcoul.kptDec.iaEnd);
#endif
            if (kcoul.kk[ik] == 0.0) continue;
            kx   = kcoul.kx[ik];
            ky   = kcoul.ky[ik];
            kz   = kcoul.kz[ik];
            //Ak   = kcoul.Ak[ik];
            for (ia=0;ia<nion;ia++){
                rxi = rx[ia];
                ryi = ry[ia];
                rzi = rz[ia];
                k_dot_r  = ( kx * rxi + ky * ryi + kz * rzi );
                kcoul.ckria[ik][ia]=cos(k_dot_r);
                kcoul.skria[ik][ia]=sin(k_dot_r);
/*
                if ( lqchtask ) {
                    qi  = q[ia];
                    rhonk_R  += qi * ckria;
                    rhonk_I  += qi * skria;
                }
                if ( ldiptask ) {
                    k_dot_mu = ( mu[ia][0] * kx + mu[ia][1] * ky + mu[ia][2] * kz );
                    rhonk_R += -k_dot_mu * skria;
                    rhonk_I +=  k_dot_mu * ckria;
                }

            } 
            kcoul.str[ik] = (rhonk_R*rhonk_R + rhonk_I*rhonk_I) * Ak;
            kcoul.rhon_R[ik] = rhonk_R;
            kcoul.rhon_I[ik] = rhonk_I;
*/
            } /* ia */ 
        } /* kpoint sum */
    } /*omp*/
}
void charge_density_q(int ik, double *q){

    double kx,ky,kz;
    double rhonk_R,rhonk_I;
    double qi;

    kx   = kcoul.kx[ik];
    ky   = kcoul.ky[ik];
    kz   = kcoul.kz[ik];
    rhonk_R = 0.0;
    rhonk_I = 0.0;
    for (int ia=0;ia<nion;ia++){
        qi=q[ia];
        rhonk_R += qi * kcoul.ckria[ik][ia];
        rhonk_I += qi * kcoul.skria[ik][ia];
    }
    kcoul.rhon_R[ik] += rhonk_R;
    kcoul.rhon_I[ik] += rhonk_I;

}

void charge_density_mu(int ik, double (*mu)[3]){

    double kx,ky,kz;
    double k_dot_mu;
    double rhonk_R,rhonk_I;

    kx   = kcoul.kx[ik];
    ky   = kcoul.ky[ik];
    kz   = kcoul.kz[ik];
    rhonk_R = 0.0;
    rhonk_I = 0.0;
    for (int ia=0;ia<nion;ia++){
        k_dot_mu = ( mu[ia][0] * kx + mu[ia][1] * ky + mu[ia][2] * kz );
        rhonk_R += -k_dot_mu * kcoul.skria[ik][ia];
        rhonk_I +=  k_dot_mu * kcoul.ckria[ik][ia];
    }
    kcoul.rhon_R[ik] += rhonk_R;
    kcoul.rhon_I[ik] += rhonk_I;

}

void charge_density_qmu(int ik, double *q, double (*mu)[3] ){

    double kx,ky,kz;
    double rhonk_R,rhonk_I;
    double qi;
    double k_dot_mu;

    kx   = kcoul.kx[ik];
    ky   = kcoul.ky[ik];
    kz   = kcoul.kz[ik];
    rhonk_R = 0.0;
    rhonk_I = 0.0;
    for (int ia=0;ia<nion;ia++){
        qi=q[ia];
        rhonk_R += qi * kcoul.ckria[ik][ia];
        rhonk_I += qi * kcoul.skria[ik][ia];
        k_dot_mu = ( mu[ia][0] * kx + mu[ia][1] * ky + mu[ia][2] * kz );
        rhonk_R += -k_dot_mu * kcoul.skria[ik][ia];
        rhonk_I +=  k_dot_mu * kcoul.ckria[ik][ia];
    }
    kcoul.rhon_R[ik] += rhonk_R;
    kcoul.rhon_I[ik] += rhonk_I;

}
