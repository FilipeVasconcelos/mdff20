#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "constants.h"
#include "config.h"
#include "field.h"
#include "integration.h"
#include "kinetic.h"
#include "md.h"
#include "io.h"
#include "tools.h"
#include "global.h"
#include "thermo.h"
#include "timing.h"
#include "verlet.h"

//#define DEBUG_MD

/******************************************************************************/
int read_md (char * controlfn) 
{
   char buffer[MAX_LEN+1];
   FILE * fp;

   fp = fopen (controlfn, "r");

   if (NULL == fp ) {
       perror("opening database file");
       return (-1);
   }

   while (EOF != fscanf(fp, "%s\n", buffer)) {
        if (strcmp(buffer,"integrator") == 0 ) {
            fscanf(fp,"%s",buffer);
            egrator=check_string("integrator",buffer,allwd_integrator_str,ALLWD_INTEGRATOR_STR); 
        } 
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
        if (strcmp(buffer,"fprint") == 0 ) {
            fscanf(fp,"%d",&fprint);
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

/******************************************************************************/
void init_md(char * controlfn){
    /* gen allowed input strings for integrator */
    strcpy(allwd_integrator_str[0 ],"nve-lf");
    strcpy(allwd_integrator_str[1 ],"nve-vv");

    /* default values */
    egrator=1; /* nve-vv */
    read_md(controlfn);
    /* tauTberendsen == dt => simple velocity rescale */ 
    if ( (nequil > 0 ) && (tauTberendsen==0.0)) tauTberendsen=dt;

    temp           = temp           * boltz_unit ; // temp = kB * T 
    dt             = dt             * time_unit  ; // angstrom*(atomicmassunit/eV)** 0.5  <= ps
    tauTberendsen  = tauTberendsen  * time_unit  ;
    if (egrator == 0 && nequil >0 ){
        lleapequi=true;
        egrator=1; /* switch to nve-vv */ 
    }

    info_md();

}

/******************************************************************************/
void info_md(){
    if (ionode) {
        SEPARATOR;
        printf("md info \n");
        LSEPARATOR;
        printf("temperature           = %.5f (K) \n", temp/boltz_unit);
        printf("dt                    = %.5f (ps)\n", dt/time_unit);
        printf("tauTberendsen         = %.5f (ps)\n", tauTberendsen/time_unit);
        printf("npas                  = %d \n", npas );
        printf("nprint                = %d \n", nprint );
        printf("nequil                = %d \n", nequil );
        putchar('\n');
    }
}


/******************************************************************************/
void run_md()
{

    calc_temp(&temp_r,&e_kin,0);

    io_node printf("in run_md\n");
    io_node printf("after previous step (leap-frog)\n"); 
    engforce();
    if(ionode){
        SEPARATOR;
        printf("properties at t=0\n");
        SEPARATOR;
        putchar('\n');
    }

    info_thermo(0,NULL); /* at t=0 */ 

#ifdef DEBUG_MD
    sample_config(0);
#endif 

    if( (ionode) && (npas>0)){
        SEPARATOR;
        printf("starting main MD loop\n");
        SEPARATOR;
        putchar('\n');
    }

    FILE *fpOSZIFF;
    fpOSZIFF = fopen ("OSZIFF", "w");
    if (NULL == fpOSZIFF ) {
        perror("opening database file");
        exit(-1);
    }

    /* ----------------------------------------------------*/
    /*                  MAIN LOOP                          */
    /* ----------------------------------------------------*/
    statime(0); /* md */
    for(istep=1; istep < npas+1; istep++) {

#ifdef DEBUG_MD
    sample_config(0);
#endif 
        /* --------------------------------- */
        /*     integration / propagators     */
        /* --------------------------------- */
        prop_agate();
        /* --------------------------------- */
        /*     check_verlet_list             */
        /* --------------------------------- */
        if (lverletL) check_verletlist();
        /* --------------------------------- */
        /*     rescale velocities in "NVE"   */
        /* --------------------------------- */
        if (istep < nequil) rescale_velocities(1);
        /* --------------------------------- */
        /*          stdout info              */
        /* --------------------------------- */
        if ( iopnode(istep,npas,nprint) ) {
            // md step timing info
            statime(1);
            mestime(&mdstepCPUtime,1,0);
            writime("MD",istep,1,0);
            statime(0);
            // md step thermodynamic info
            info_thermo(0,NULL); /*STDOUT*/
        }
        /* --------------------------------- */
        /*          OSZIFF info              */
        /* --------------------------------- */
        if ( iopnode(istep,npas,fprint) ) {
            info_thermo(1,fpOSZIFF);   /*OSZIFF*/
        }
    }
    /* ----------------------------------------------------*/
    write_config();
    fclose(fpOSZIFF);
}

