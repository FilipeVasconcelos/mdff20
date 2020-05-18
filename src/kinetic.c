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

void init_velocities()
{
    maxwellboltzmann_velocities();
}

void print_velocities()
{
    if (ionode){
        for ( int ia=0; ia < nion; ia++)
        {
            printf("%5d %5.3f %5.3f %5.3f\n",ia,vx[ia],vy[ia],vz[ia]);
        }
    }
}

void maxwellboltzmann_velocities()
{
    double rtemp;
    double sx=0.0; double sy=0.0; double sz=0.0;

    for ( int ia=0; ia < nion; ia++ )
    {
        rtemp = sqrt((double)temp/massia[ia]);
        //printf("rtemp %12.8f\n",rtemp);
        vx [ia] = rtemp * box_muller(0.0,1.0);
        vy [ia] = rtemp * box_muller(0.0,1.0);
        vz [ia] = rtemp * box_muller(0.0,1.0);
        //printf("boxmuller %12.8f\n",box_muller(0.0,1.0));
        sx = sx + vx[ia];
        sy = sy + vy[ia];
        sz = sz + vz[ia];
    }
    sx*=onenion;
    sy*=onenion;
    sz*=onenion;
    for ( int ia=0; ia < nion; ia++ )
    {   
        vx [ia] -= sx;
        vy [ia] -= sy;
        vz [ia] -= sz;
    }
    double temp,ekin;
    calc_temp(&temp, &ekin,0);
    io_node printf("(after maxwellboltzmann) temp : %f kin : %f \n",temp,ekin);
    temp_r= temp;
    e_kin = ekin  ;
}

void calc_temp(double *temp, double *ekin, int flag)
{
    double l = 3.0 * ((double) nion);
    double kin;
    kin = 0.0;
    for (int ia=0; ia< nion; ia++)
    {
        kin += ( vx[ia]*vx[ia] + vy[ia]*vy[ia] + vz[ia]*vz[ia] ) * massia[ia] ;
    }
    (*ekin)=(kin*0.5);
    (*temp)= (2.0 * kin) / ( l * boltz_unit );
}

void rescale_velocities()
{
    double sx=0.0; double sy=0.0; double sz=0.0;
    double T,ekin;
    calc_temp(&T,&ekin,0);
    ekin*=0.5 ;
    double lambda = sqrt(1.0 + (dt / tauTberendsen) * (( temp / T / boltz_unit ) - 1.0 )) ;
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

    calc_temp(&T,&ekin,0);
    //if (ionode && istep%nprint==0) printf("(after rescaling)        temp : %f kin : %f lambda %f\n",T,ekin,lambda);
}
