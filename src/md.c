#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"
#include "config.h"
#include "field.h"
#include "integration.h"
#include "kinetic.h"
#include "md.h"
#include "thermo.h"

int read_md (char * controlfn) 
{
   char buffer[15];
   double data;
   FILE * fp;

   fp = fopen (controlfn, "r");

   if (NULL == fp ) {
       perror("opening database file");
       return (-1);
   }

   while (EOF != fscanf(fp, "%s %lf\n", buffer,&data)) {
        if (strcmp(buffer,"dt") == 0 ) {
            dt=data;
        } 
        if (strcmp(buffer,"temp") == 0 ) {
            temp=data;
        } 
        if (strcmp(buffer,"npas") == 0 ) {
            npas=(int) data;
        } 
        if (strcmp(buffer,"nprint") == 0 ) {
            nprint=(int) data;
        } 
        if (strcmp(buffer,"nequil") == 0 ) {
            nequil=(int) data;
        } 
        if (strcmp(buffer,"tauTberendsen") == 0 ) {
            tauTberendsen=data;
        } 
   }


   fclose(fp);
   
   return(0);
}

void init_md(char * controlfn){

    read_md(controlfn);
    temp           = temp           * boltz_unit ; // temp = kB * T 
    dt             = dt             * time_unit  ; // angstrom*(atomicmassunit/eV)** 0.5  <= ps
    tauTberendsen  = tauTberendsen  * time_unit  ;
    info_md();

}

void info_md(){
    printf("temperature     : %.5f (K) \n", temp);
    printf("dt              : %.5f (ps)\n", dt );
    printf("tauTberendsen   : %.5f (ps)\n", tauTberendsen );
    printf("npas            : %d \n", npas );
    printf("nprint          : %d \n", nprint );
    printf("nequil          : %d \n", nequil );

}


void run_md()
{
    
    /*
    for(int ia=0; ia<nion; ia++) {
        printf("(in run_md) reading %d %s %lf %lf %lf\n",ia,atype[ia],rx[ia],ry[ia],rz[ia]);
    }
*/
    /*
    t=0
    */
    /*
    for(int ia=0;ia<nion;ia++) {
            rxs [ia]  = rx [ia] - vx[ia] * dt;
            rys [ia]  = ry [ia] - vy[ia] * dt;
            rzs [ia]  = rz [ia] - vz[ia] * dt;
    }
    */
    engforce();
    write_thermo();

    for(itime=1; itime<npas+1;itime++)
    {
        prop_velocity_verlet();

        if (itime < nequil) rescale_velocities();
        if (itime%nprint==0 || itime == npas ) write_thermo();
    }
}

