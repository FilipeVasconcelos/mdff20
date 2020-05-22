#include<stdlib.h>
#include<stdio.h>
#include "config.h"
#include "md.h"
#include "kinetic.h"
#include "field.h"
#include "thermo.h"
#include "timing.h"

void prop_velocity_verlet(){

    statime(12);
    double *fxs, *fys, *fzs;
    fxs=malloc(nion*sizeof(*fxs));
    fys=malloc(nion*sizeof(*fys));
    fzs=malloc(nion*sizeof(*fzs));
    double dtsq2 = dt * dt * 0.5L; 
    double dt2 = dt * 0.5L;

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
    for(int ia=0;ia<nion;ia++){
        vx [ia] += ( fxs[ia] + fx[ia] ) * dt2 * invemassia[ia];
        vy [ia] += ( fys[ia] + fy[ia] ) * dt2 * invemassia[ia];
        vz [ia] += ( fzs[ia] + fz[ia] ) * dt2 * invemassia[ia];
    }
    // -------------------------
    // full t+dt kinetic energy 
    // -------------------------
    double tempi,kin;
    calc_temp(&tempi, &kin,0) ;
    temp_r= tempi             ;
    e_kin = kin               ;
    free(fxs);
    free(fys);
    free(fzs);
    statime(13);
    mestime(&propagatorCPUtime,13,12);
}


void prop_leap_frog(){
    
    engforce();

    statime(12);
    double dtsq = dt * dt ;
    double *urx, *ury, *urz;
    double idt = 0.5 / dt ;

    urx=malloc(nion*sizeof(*urx));
    ury=malloc(nion*sizeof(*ury));
    urz=malloc(nion*sizeof(*urz));
   
    for (int ia=0; ia < nion ; ia++)    { 
        // r(t+dt) = 2 r(t) - r (t-dt) + f(t)/m dt*dt 
        urx [ia] = 2.0 * rx [ ia ] - rxs [ ia ] + fx [ ia ] * dtsq * invemassia[ia];
        ury [ia] = 2.0 * ry [ ia ] - rys [ ia ] + fy [ ia ] * dtsq * invemassia[ia];
        urz [ia] = 2.0 * rz [ ia ] - rzs [ ia ] + fz [ ia ] * dtsq * invemassia[ia];

        //v(t) = ( r(t+dt) - r(t) ) / ( 2 dt)
        vx [ia] = idt * ( urx [ ia ] - rxs [ ia ] )  ;
        vy [ia] = idt * ( ury [ ia ] - rys [ ia ] )  ;
        vz [ia] = idt * ( urz [ ia ] - rzs [ ia ] )  ;
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
    double tempi,kin;
    calc_temp ( &tempi, &kin,0) ;
    temp_r= tempi             ;
    e_kin = kin               ;
    free(urx);
    free(ury);
    free(urz);
    statime(13);
    mestime(&propagatorCPUtime,13,12);

}

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

    switch (egrator)
    {
        case 0: /* nve-lf */
            prop_leap_frog();
            break;
        case 1: /* nve-vv */
            prop_velocity_verlet();
            break;
    }
}


