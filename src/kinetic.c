#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "rand.h"
#include "config.h"
#include "md.h"
#include "kinetic.h"

void init_velocities()
{
    printf("inside init_velocities %d\n",nion);
    for ( int ia=0; ia < nion; ia++)
    {
        vx[ia]=0.0;
        vy[ia]=0.0;
        vz[ia]=0.0;
    }

}
void print_velocities()
{
    printf("inside print_velocities %d\n",nion);
    for ( int ia=0; ia < nion; ia++)
    {
        printf("%5d %5.3f %5.3f %5.3f\n",ia,vx[ia],vy[ia],vz[ia]);
    }
    printf("out print_velocities %d\n",nion);

}

void maxwellboltzmann_velocities()
{
    printf("inside maxwellboltzmann_velocities %d\n",nion);
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
    printf("(after maxwell) temp : %f kin : %f \n",temp,ekin);
    printf("here\n");

}

void calc_temp(double *temp, double *ekin, int flag)
{
    double l = 3.0 * (double) nion;
    double kin;

    kin = 0.0;
    for (int ia=0; ia< nion; ia++)
    {
        kin =  kin + ( vx[ia]*vx[ia] + vy[ia]*vy[ia] + vz[ia]*vz[ia] ) * massia[ia] ;
    //    printf("(in calc_temp %d %f %f %f %f %f\n",ia,kin, vx[ia], vy[ia], vz[ia],massia[ia]);
    }
    *ekin=(kin*0.5);
    *temp = (2.0 * kin / ( l * boltz_unit ));
//    printf("out calc_temp\n");
//    if (flag) exit(0);
    

}

void rescale_velocities()
{
    double sx=0.0; double sy=0.0; double sz=0.0;
    double T,ekin;
    
    calc_temp(&T,&ekin,0);
    //printf("(before rescaling) temp : %f kin : %f \n",T,ekin);
   
    double lambda = sqrt( ( 1.0 + (dt / tauTberendsen) * ( ( temp / T / boltz_unit ) - 1.0 ) ) ) ;
    //printf("lambda %f dt %f tau %f temp %f T %f boltz %f\n",lambda,dt,tauTberendsen,temp/boltz_unit,T/boltz_unit,boltz_unit); 

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
    if (itime%nprint==0) printf("(after rescaling) temp : %f kin : %f lambda %f\n",T,ekin,lambda);
}




























