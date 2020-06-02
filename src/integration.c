#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "config.h"
#include "md.h"
#include "io.h"
#include "kinetic.h"
#include "field.h"
#include "thermo.h"
#include "timing.h"
#include "integration.h"

/******************************************************************************/
void prop_velocity_verlet(){

    statime(12);
    double *fxs, *fys, *fzs;
    fxs=malloc(nion*sizeof(*fxs));
    fys=malloc(nion*sizeof(*fys));
    fzs=malloc(nion*sizeof(*fzs));
    double dtsq2 = dt * dt * 0.5; 
    double dt2 = dt * 0.5;

    for(int ia=0;ia<nion;ia++){
        fxs[ia]=fx[ia];fys[ia]=fy[ia];fzs[ia]=fz[ia];
        rx[ia] += vx[ia] * dt + (fx[ia] * dtsq2 ) * invemassia[ia] ;
        ry[ia] += vy[ia] * dt + (fy[ia] * dtsq2 ) * invemassia[ia] ;
        rz[ia] += vz[ia] * dt + (fz[ia] * dtsq2 ) * invemassia[ia] ;
    }
    statime(13);
    mestime(&propagatorCPUtime,13,12);
// ------------------
    engforce();
// ------------------
    statime(12);
    double kin=0.0;
    for(int ia=0;ia<nion;ia++){
        vx [ia] += ( fxs[ia] + fx[ia] ) * dt2 * invemassia[ia];
        vy [ia] += ( fys[ia] + fy[ia] ) * dt2 * invemassia[ia];
        vz [ia] += ( fzs[ia] + fz[ia] ) * dt2 * invemassia[ia];
        kin     += (vx [ia]* vx [ia]+ vy [ia]*vy [ia] + vz [ia]*vz [ia] )*massia[ia];
    }
    // -------------------------
    // full t+dt kinetic energy 
    // -------------------------
    e_kin=kin*0.5;
    free(fxs);
    free(fys);
    free(fzs);
    statime(13);
    mestime(&propagatorCPUtime,13,12);
}

/******************************************************************************/
void prop_pos_vel_verlet(double* kin){

    double dt2 = dt * 0.5;

    for(int ia=0;ia<nion;ia++){
        vx [ia] += fx[ia] * dt2 * invemassia[ia];
        vy [ia] += fy[ia] * dt2 * invemassia[ia];
        vz [ia] += fz[ia] * dt2 * invemassia[ia];
        rx[ia] += vx[ia] * dt ;
        ry[ia] += vy[ia] * dt ;
        rz[ia] += vz[ia] * dt ;
    }

    /********************/
    engforce();
    /********************/
    *kin=0.0;
    for(int ia=0;ia<nion;ia++){
        vx [ia] += fx[ia]  * dt2 * invemassia[ia];
        vy [ia] += fy[ia]  * dt2 * invemassia[ia];
        vz [ia] += fz[ia]  * dt2 * invemassia[ia];
        *kin    += (vx[ia] * vx[ia] + vy[ia] * vy[ia] + vz[ia] * vz[ia] ) * massia[ia];
    }
    *kin*=0.5;
}

/******************************************************************************/
void prop_leap_frog(){
    
    engforce();

    statime(12);
    double dtsq = dt * dt ;
    double *urx, *ury, *urz;
    double idt = 0.5 / dt ;

    urx=malloc(nion*sizeof(*urx));
    ury=malloc(nion*sizeof(*ury));
    urz=malloc(nion*sizeof(*urz));
  
    double kin=0.0; 
    for (int ia=0; ia < nion ; ia++)    { 
        // r(t+dt) = 2 r(t) - r (t-dt) + f(t)/m dt*dt 
        urx [ia] = 2.0 * rx [ ia ] - rxs [ ia ] + fx [ ia ] * dtsq * invemassia[ia];
        ury [ia] = 2.0 * ry [ ia ] - rys [ ia ] + fy [ ia ] * dtsq * invemassia[ia];
        urz [ia] = 2.0 * rz [ ia ] - rzs [ ia ] + fz [ ia ] * dtsq * invemassia[ia];

        //v(t) = ( r(t+dt) - r(t) ) / ( 2 dt)
        vx [ia] = idt * ( urx [ ia ] - rxs [ ia ] )  ;
        vy [ia] = idt * ( ury [ ia ] - rys [ ia ] )  ;
        vz [ia] = idt * ( urz [ ia ] - rzs [ ia ] )  ;
        kin += (vx [ia]*vx [ia] + vy [ia] * vy [ia] + vz [ia]*vz [ia]) *massia[ia];
        // updated positions r(t-dt) <= r(t)  and r(t) <= r (t+dt) 
        rxs [ia] = rx  [ia] ;  
        rys [ia] = ry  [ia] ;  
        rzs [ia] = rz  [ia] ;  
        rx  [ia] = urx [ia] ;  
        ry  [ia] = ury [ia] ;  
        rz  [ia] = urz [ia] ;  
    }
    // -------------------------
    // full t+dt kinetic energy 
    // -------------------------
    kin*=0.5;
    e_kin=kin;
    free(urx);
    free(ury);
    free(urz);
    statime(13);
    mestime(&propagatorCPUtime,13,12);

}

/******************************************************************************/
void nhcn(){

    double L = 3.0 * (double) nion ;
    double *Q;
    double kini;

    Q=malloc(nhc_n*sizeof(*Q));
    for(int i=0;i<nhc_n;i++){
        Q[i]= timesca_thermo*timesca_thermo* temp;
    }
    Q[0]=Q[0]*L;

    kini=calc_kin();
    chain_nhcn(&kini, vxi, xi, Q, L);
    prop_pos_vel_verlet(&kini);
    chain_nhcn(&kini, vxi, xi, Q, L);
    e_kin=calc_kin();
    //e_kin=kini*0.5;

    // ===========================================
    //  conserved quantity energy in NVT ensemble
    // ===========================================
    // note on units :
    // [ vxi ] = [ eV ] [ T ]
    // [ Q ]   = [ eV ] [ T ]**2
    // vxi * vxi / Q =  [ eV ]
    // [ temp ] * [xi] = [ eV ]  ( [xi] = sans unité )
    double nvt1,nvt2,nvt3,nvt4;
    nvt1 = L * temp * xi[0];
    nvt2 = vxi[0] * vxi[0] * 0.5 / Q[0];
    nvt3 = 0.0;
    nvt4 = 0.0;
    for (int inhc = 1;inhc<nhc_n;inhc++){
	nvt3 += vxi[inhc] * vxi[inhc] * 0.5 / Q[inhc];
        nvt4 += temp * xi[inhc];
    }
    e_nvt = nvt1 + nvt2 + nvt3 + nvt4;
    io_pnode printf("  e_nvt"ff5"\n",e_nvt,nvt1,nvt2,nvt3,nvt4);

    free(Q);
}

/******************************************************************************/
// [1] Molecular Physics, (1996), v87, n5, p1117 Martyna and al.
// [2] Phys Lett. A, (1190), v150 n5,6,7, p262, Yoshida
// [3] https://files.nyu.edu/mt33/public/abstracts/a6_19_s18.pdf
/******************************************************************************/
void chain_nhcn ( double *kin , double vxi[nhc_n] , double xi[nhc_n] , double Q[nhc_n] , double L ){

    //printf("in chain_nhcn\n");
    double *yosh_w, *dt_yosh, *G;
    yosh_w = malloc(nhc_yosh_order *sizeof(*yosh_w));
    dt_yosh = malloc(nhc_yosh_order *sizeof(*dt_yosh));
    G = malloc(nhc_n*sizeof(*G));
    for (int i=0;i<nhc_n;i++){
        G[i]=0.0;
    }
    for (int i=0;i<nhc_yosh_order;i++){
        yosh_w[i]=0.0;
        dt_yosh[i]=0.0;
    }
    //printf("here!!!\n");

    switch ( nhc_yosh_order ) {
        default : 
            io_node printf("value of yoshida order not available");
            break;
        case 1:
            yosh_w[0]=1.0;
            break;
        case 3:
            yosh_w[0] = 1.0/(2.0-pow(2.0,1.0/3.0));
            yosh_w[1] = 1.0 - 2.0 * yosh_w[0] ;
            yosh_w[2] = yosh_w[0];
            break;
        case 5:
            yosh_w[0] = 1.0/(4.0-pow(4,1.0/3.0));
            yosh_w[1] = yosh_w[0];
            yosh_w[2] = yosh_w[0];
            yosh_w[3] = yosh_w[0];
            yosh_w[4] = 1.0 - 4.0 * yosh_w[0];
            break;
        case 7:
            yosh_w[0]=  .78451361047756;
            yosh_w[1]=  .235573213359357;
            yosh_w[2]=-1.17767998417887;
            yosh_w[3]= 1.0-2.0*(yosh_w[0]+yosh_w[1]+yosh_w[2]);
            yosh_w[4]=yosh_w[2];
            yosh_w[5]=yosh_w[1];
            yosh_w[6]=yosh_w[0];
            break;
        case 9:
            break;
    }
    for (int i=0;i<nhc_yosh_order;i++){
        dt_yosh[i]=yosh_w[i]*dt/((double) nhc_mults);
    }
    double s=1.0;
    G[0] =  2.0 * (*kin) - L * temp;

    double dts,dts2,dts4,dts8;

    for (int k=0;k<nhc_mults;k++){
        for (int j=0;j<nhc_yosh_order;j++){
            dts  = dt_yosh[j];
            dts2 = dts  * 0.5;
            dts4 = dts2 * 0.5;
            dts8 = dts4 * 0.5;
            G[nhc_n-1] = ( vxi[nhc_n-2] * vxi[nhc_n-2] / Q[nhc_n-2] - temp);
            vxi[nhc_n-1] += G[nhc_n-1] * dts4;
            for (int inh=nhc_n-2;inh>-1;inh--){
                //scale thermo momentum
                vxi[inh] *= exp ( -vxi [inh+1] * dts8 / Q [inh+1] );
                //propagate thermo momentum
                vxi[inh] += G[inh]*dts4;
                //scale thermo momentum
                vxi[inh] *= exp ( -vxi [inh+1] * dts8 / Q [inh+1] );
            }
            //exp5: scale velocities
            s *= exp ( - vxi[0] * dts2 / Q [0] );
            //exp6 : propagating xi
            //minus in [3] seems to be wrong ????
            for(int i=0;i<nhc_n;i++){
                xi[i] += vxi[i] * dts2 / Q[i];
            }
            G[0] = 2.0*(*kin)*s*s - L * temp;
            for (int inh=0;inh<nhc_n-2;inh++){
                //scale thermo momentum
                vxi[inh] *= exp ( -vxi [inh+1] * dts8 / Q [inh+1] );
                //propagate thermo momentum
                vxi[inh] += G[inh]*dts4;
                //scale thermo momentum
                vxi[inh] *= exp ( -vxi [inh+1] * dts8 / Q [inh+1] );
                G[inh+1] = ( vxi [inh] * vxi[inh] / Q[inh] - temp);
            }
            vxi[nhc_n-1] += G[nhc_n-1] * dts4;
        }
    }
    *kin=0.0;
    for(int ia=0; ia<nion; ia++) {
        vx [ia] *= s;
        vy [ia] *= s;
        vz [ia] *= s;
        *kin    += ( vx[ia] * vx [ia] + vy[ia] * vy[ia] + vz[ia] * vz[ia] ) * massia[ia];
    }

    *kin*=0.5;
    free(G);
    free(yosh_w);
    free(dt_yosh);

}


/******************************************************************************/
void prop_agate(int step){

    statime(12);
    /* ==========================================================
       leap-frog algorithm set the first leap value for r(t-dt) 
       if first steps are equilibrated with a velocity-verlet
      ==========================================================*/
/*    if ( lleapequi && istep >= nequil ) {
        egrator = 0;
        lleapequi = false;
        for(int ia=0;ia<nion;ia++) {
            rxs [ia]  = rx [ia] - vx[ia] * dt;
            rys [ia]  = ry [ia] - vy[ia] * dt;
            rzs [ia]  = rz [ia] - vz[ia] * dt;
        }
        printf("  switch to leap-frog %d\n",egrator);    
    }
    */
    statime(13);
    mestime(&propagatorCPUtime,13,12);
    switch (egrator) {
        case 0: /* nve-lf */
            prop_leap_frog();
            break;
        case 1: /* nve-vv */
            prop_velocity_verlet();
            break;
        case 2: /* nvt-nhc_n*/
            nhcn();
            break;
    }

}


