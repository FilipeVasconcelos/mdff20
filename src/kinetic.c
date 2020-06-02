#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "io.h"
#include "rand.h"
#include "config.h"
#include "thermo.h"
#include "kinetic.h"
#include "md.h"
#include "tools.h"

/******************************************************************************/
void init_velocities()
{
    double comit[NTYPEMAX+1][3];
    

    maxwellboltzmann_velocities();
    rescale_velocities(0);

    com(vx,vy,vz,ntype+1, comit);
    if (ionode) {
        //printf("c.o.m. vel ALL "ee3"\n",comit[ntype][0],comit[ntype][1],comit[ntype][2]);
        for (int it=0;it<ntype;it++){
        //    printf("c.o.m. vel %3s "ee3"\n",atypit[it],comit[it][0],comit[it][1],comit[it][2]);
        }
    }
    sample_config(0);
//    exit(-1);
}

/******************************************************************************/
void print_velocities()
{
    if (ionode){
        for ( int ia=0; ia < nion; ia++)
        {
            printf("%5d %5.3f %5.3f %5.3f\n",ia,vx[ia],vy[ia],vz[ia]);
        }
    }
}

/******************************************************************************/
void maxwellboltzmann_velocities()
{
    double rtemp;
    double sx=0.0; double sy=0.0; double sz=0.0;

    for ( int ia=0; ia < nion; ia++ ) {
        rtemp = sqrt(temp*invemassia[ia]);
//        printf("rtemp %12.8f\n",rtemp);
        vx [ia] = rtemp * box_muller(0.0,1.0);
        vy [ia] = rtemp * box_muller(0.0,1.0);
        vz [ia] = rtemp * box_muller(0.0,1.0);
        //printf("boxmuller %12.8f\n",box_muller(0.0,1.0));
        sx += vx[ia];
        sy += vy[ia];
        sz += vz[ia];
    }
    sx*=onenion;
    sy*=onenion;
    sz*=onenion;
    for ( int ia=0; ia < nion; ia++ ) {   
        vx [ia] += -sx;
        vy [ia] += -sy;
        vz [ia] += -sz;
    }
}


/******************************************************************************/
double calc_kin() {
    double kin = 0.0;
    for (int ia=0; ia< nion; ia++) {
        kin += ( vx[ia]*vx[ia] + vy[ia]*vy[ia] + vz[ia]*vz[ia] ) * massia[ia] ;
    }
    kin*=0.5;
    return kin;
}

/******************************************************************************/
double calc_temp(double kin) {
    return ( 2.0 * kin ) / ( 3.0 * ((double) nion)  * boltz_unit );
}

/******************************************************************************/
void rescale_velocities(int quiet)
{
    double sx=0.0; double sy=0.0; double sz=0.0;
    double l = 3.0 * ((double) nion);
    double ekin = calc_kin();
    double T    = calc_temp(ekin);

    double lambda;
    if (egrator !=2 ) { 
        lambda = sqrt(1.0 + (dt / tauTberendsen) * (( temp / T / boltz_unit ) - 1.0 )) ;
    } 
    else {
        lambda = sqrt(1.0 + (( temp / T / boltz_unit ) - 1.0 )) ;
    }
    for(int ia=0; ia<nion; ia++) {
        vx [ia] *= lambda;
        vy [ia] *= lambda;
        vz [ia] *= lambda;
        sx += vx [ia] ;
        sy += vy [ia] ;
        sz += vz [ia] ;
    }
    sx*=onenion;
    sy*=onenion;
    sz*=onenion;
    for(int ia=0; ia<nion; ia++) {
        vx [ia] -= sx;
        vy [ia] -= sy;
        vz [ia] -= sz;
    }
    if (ionode && quiet == 0) {
        printf("effective temperature T  = %10.4f\n",T);
        printf("wanted temperature    T0 = %10.4f\n",temp/boltz_unit);
        if (egrator !=2 ) {
            printf("velocities rescaled by   : "ee ee"\n",lambda,(dt/tauTberendsen)*((temp/T/boltz_unit)-1.0));
        }
        else {
            printf("velocities rescaled by   : "ee ee"\n",lambda,(temp/T/boltz_unit)-1.0);
        }
    }
    ekin = calc_kin();
    T    = calc_temp(ekin);
    if (ionode && istep%nprint==0) {
        printf("(after rescaling)\n");
        printf("temp                     : %f\n",T);
        printf("kin                      : %f\n",ekin);
        printf("lambda                   : %f\n",lambda);
    }
}
