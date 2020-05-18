#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"
#include "config.h"
#include "field.h"
#include "integration.h"
#include "kinetic.h"
#include "md.h"
#include "io.h"
#include "thermo.h"
#include "timing.h"

int read_md (char * controlfn) 
{
   char buffer[15];
   FILE * fp;

   fp = fopen (controlfn, "r");

   if (NULL == fp ) {
       perror("opening database file");
       return (-1);
   }

   while (EOF != fscanf(fp, "%s\n", buffer)) {
        if (strcmp(buffer,"dt") == 0 ) {
            fscanf(fp,"%lf",&dt);
        } 
        if (strcmp(buffer,"temp") == 0 ) {
            fscanf(fp,"%lf",&temp);
        } 
        if (strcmp(buffer,"npas") == 0 ) {
            fscanf(fp,"%d",&npas);
        } 
        if (strcmp(buffer,"nprint") == 0 ) {
            fscanf(fp,"%d",&nprint);
        } 
        if (strcmp(buffer,"nequil") == 0 ) {
            fscanf(fp,"%d",&nequil);
        } 
        if (strcmp(buffer,"tauTberendsen") == 0 ) {
            fscanf(fp,"%lf",&tauTberendsen);
        } 
   }
   fclose(fp);
   return 0;
}

void init_md(char * controlfn){

    read_md(controlfn);
    info_md();
    temp           = temp           * boltz_unit ; // temp = kB * T 
    dt             = dt             * time_unit  ; // angstrom*(atomicmassunit/eV)** 0.5  <= ps
    tauTberendsen  = tauTberendsen  * time_unit  ;

}

void info_md(){
    if (ionode) {
        SEPARATOR;
        printf("md info \n");
        LSEPARATOR;
        printf("temperature           = %.5f (K) \n", temp);
        printf("dt                    = %.5f (ps)\n", dt);
        printf("tauTberendsen         = %.5f (ps)\n", tauTberendsen);
        printf("npas                  = %d \n", npas );
        printf("nprint                = %d \n", nprint );
        printf("nequil                = %d \n", nequil );
        putchar('\n');
    }
}


void run_md()
{
    /* previous step in leap-frog integration */ 
    for(int ia=0;ia<nion;ia++) {
            rxs [ia]  = rx [ia] - vx[ia] * dt;
            rys [ia]  = ry [ia] - vy[ia] * dt;
            rzs [ia]  = rz [ia] - vz[ia] * dt;
    }
    
    engforce();
    if(ionode){
        SEPARATOR;
        printf("properties at t=0\n");
        SEPARATOR;
        putchar('\n');
    }
    info_thermo(); /* at t=0 */ 
    if( (ionode) && (npas>0)){
        SEPARATOR;
        printf("starting main MD loop\n");
        SEPARATOR;
        putchar('\n');
    }

    /* ----------------------------------------------------*/
    /*                  MAIN LOOP                          */
    /* ----------------------------------------------------*/
    statime(0);
    for(istep=1; istep < npas+1; istep++) {
        
        // integration / propagators
        //prop_leap_frog();
        prop_velocity_verlet();
        
        if (istep < nequil) rescale_velocities();
        
        if (istep % nprint==0 || istep == npas ) {
            statime(1);
            info_thermo();
            mestime(&mdstepCPUtime,1,0);
            writime("MD",istep,1,0);
            statime(0);
        }

    }
    /* ----------------------------------------------------*/
}

