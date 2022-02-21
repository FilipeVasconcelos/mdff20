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

#ifdef DEBUG
    #define DEBUG_MD
#endif
//#define DEBUG_
#ifdef DEBUG_
    #define DEBUG_MD
#endif

/******************************************************************************/
int read_md (char * controlfn)
{
    char buffer[MAX_LEN+1];
    FILE * fp;

    fp = fopen (controlfn, "r");

    if (NULL == fp ) {
        pError("opening control file (reading md) ");
        return (-1);
    }

    while (EOF != fscanf(fp, "%s\n", buffer)) {
        if (strcmp(buffer,"integrator") == 0 ) {
            fscanf(fp,"%s",buffer);
            egrator=check_string("integrator",buffer,allwd_integrator,ALLWD_INTEGRATOR_STR);
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
        if (strcmp(buffer,"cprint") == 0 ) {
            fscanf(fp,"%d",&cprint);
        }
        if (strcmp(buffer,"nequil") == 0 ) {
            fscanf(fp,"%d",&nequil);
        }
        if (strcmp(buffer,"nequilT") == 0 ) {
            fscanf(fp,"%d",&nequilT);
        }
        if (strcmp(buffer,"tauTberendsen") == 0 ) {
            fscanf(fp,"%lf",&tauTberendsen);
        }
        if (strcmp(buffer,"nhc_n") == 0 ) {
            fscanf(fp,"%d",&nhc_n);
        }
        if (strcmp(buffer,"nhc_mults") == 0 ) {
            fscanf(fp,"%d",&nhc_mults);
        }
        if (strcmp(buffer,"nhc_yosh_order") == 0 ) {
            fscanf(fp,"%d",&nhc_yosh_order);
        }
        if (strcmp(buffer,"timesca_thermo") == 0 ) {
            fscanf(fp,"%lf",&timesca_thermo);
        }
        if (strcmp(buffer,"timesca_baro") == 0 ) {
            fscanf(fp,"%lf",&timesca_baro);
        }
   }
   fclose(fp);
   return 0;
}
/******************************************************************************/
void default_md(){
    /* gen allowed input strings for integrator */
    strcpy(allwd_integrator[0],"nve-lf");
    strcpy(allwd_integrator[1],"nve-vv");
    strcpy(allwd_integrator[2],"nvt-nhcn");
    /* gen allowed input int for yosh param 1,3,5,7,9*/
    for (int i =0;i<ALLWD_YOSH_PARAM;i++){
        yosh_allowed[i]=2*i+1;
    }
    /* gen allowed rescale integrator */
    allwd_rescale_integrator[0]=1; /*nve-vv*/
    egrator=1; /* nve-vv */

    nprint=1;
    fprint=1;
    cprint=1000;
}
/******************************************************************************/
void check_md(){
    /* tauTberendsen == dt => simple velocity rescale */
    if ( (nequil > 0 ) && (tauTberendsen==0.0)) tauTberendsen=dt;

    temp           = temp           * boltz_unit ; // temp = kB * T
    dt             = dt             * time_unit  ; // angstrom*(atomicmassunit/eV)** 0.5  <= ps
    tauTberendsen  = tauTberendsen  * time_unit  ;
    timesca_thermo = timesca_thermo * time_unit  ;
    if (egrator == 0 && nequil >0 ){
        lleapequi=true;
        egrator=1; /* switch to nve-vv */
    }
    rescale_allowed=false;
    for (int i =0;i<ALLWD_RESCALE_INTEGRATOR;i++){
        if ( egrator == allwd_rescale_integrator[i] ) rescale_allowed=true;
    }
}

/******************************************************************************/
void init_md(char * controlfn){

    /* default values */
    default_md();
    /* read parameters */
    read_md(controlfn);
    /* check parameters */
    check_md();
    /* alloc memory */
    alloc_md();
    /* print information */
    info_md();

}

void alloc_md(){
    if (egrator == 2 ){
        vxi=malloc(nhc_n*sizeof(*vxi));
        xi=malloc(nhc_n*sizeof(*xi));
    }
}

void free_md(){
    if (egrator == 2 ){
        free(vxi);
        free(xi);
    }
}

/******************************************************************************/
void info_md(){
    if (ionode) {
        SEPARATOR;
        printf("md info \n");
        LSEPARATOR;
        printf("integrator            = %s\n", allwd_integrator[egrator]);
        printf("temperature           = %.5f (K) \n", temp/boltz_unit);
        printf("dt                    = %.5f (ps)\n", dt/time_unit);
        printf("tauTberendsen         = %.5f (ps)\n", tauTberendsen/time_unit);
        printf("npas                  = %d \n", npas );
        printf("nprint                = %d \n", nprint );
        printf("nequil                = %d \n", nequil );
        printf("nequilT               = %d \n", nequilT );
        putchar('\n');
    }
}


/******************************************************************************/
void run_md()
{

    e_kin=calc_kin();

#ifdef DEBUG_MD
    sample_config(0);
#endif
    
    if(ionode){
        SEPARATOR;
        printf("properties at t=0\n");
        SEPARATOR;
        putchar('\n');
    }
    /* energy, field, forces */
    engforce();

    info_thermo(0,NULL); /* at t=0 */

    if( (ionode) && (npas>0)){
        SEPARATOR;
        printf("starting main MD loop\n");
        SEPARATOR;
        putchar('\n');
    }
    FILE *fpOSZIFF;
    fpOSZIFF = fopen ("OSZIFF", "w");
    if (NULL == fpOSZIFF ) {
        pError("opening OSZIFF file");
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
        prop_agate(istep);
        /* --------------------------------- */
        /*     check_verlet_list             */
        /* --------------------------------- */
        if (lverletL) check_verletlist();
        /* --------------------------------- */
        /*     rescale velocities in "NVE"   */
        /* --------------------------------- */
        if ( ( rescale_allowed    ) &&
             ( istep < nequil     ) &&
             ( istep%nequilT == 0 ) )  rescale_velocities(1);
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
        /* --------------------------------- */
        /*          CONTFF info              */
        /* --------------------------------- */
        if ( iopnode(istep,npas,cprint) ) {
            write_config();
        }
    }
    /* ----------------------------------------------------*/
    write_config();
    fclose(fpOSZIFF);
}

