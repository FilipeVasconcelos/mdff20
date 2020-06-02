#include <stdio.h>
#include <string.h>

#include "constants.h"
#include "config.h"
#include "thermo.h"
#include "field.h"
#include "nmlj.h"
#include "ewald.h"
#include "cell.h"
#include "timing.h"
#include "md.h"
#include "io.h"
#include "tools.h"

/******************************************************************************/
int read_field(char* controlfn)
{
   char buffer[MAX_LEN+1];
   FILE * fp;
   fp = fopen (controlfn, "r");
   if (NULL == fp )  {
       perror("opening database file");
       return (-1);
   }
   while (EOF != fscanf(fp, "%s\n", buffer)) { 
        if (strcmp(buffer,"lnmlj") == 0 ) {
            fscanf(fp,"%s",buffer);
            lnmlj=check_boolstring("lnmlj",buffer); 
        } 
        if (strcmp(buffer,"lcoul") == 0 ) {
            fscanf(fp,"%s",buffer);
            lnmlj=check_boolstring("lcoul",buffer); 
        } 
        // mass of type
        if (strcmp(buffer,"massit") == 0 ) {
            for(int it=0;it<ntype;it++){
                fscanf(fp,"%lf",&massit[it]);
            }
        } 
        // charges on type
        if (strcmp(buffer,"qit") == 0 ) {
            for(int it=0;it<ntype;it++){
                fscanf(fp,"%lf",&qit[it]);
            }
        } 
        // dipole on type
        if (strcmp(buffer,"dipit") == 0 ) {
            for(int it=0;it<ntype;it++){
                fscanf(fp,"%lf",&dipit[it]);
            }
        } 
        // quadrupole on type
        if (strcmp(buffer,"quadit") == 0 ) {
            for(int it=0;it<ntype;it++){
                fscanf(fp,"%lf",&quadit[it]);
            }
        } 
        if (strcmp(buffer,"alphaES") == 0 ) {
            fscanf(fp,"%lf",&alphaES);
        } 

   }
   fclose(fp);
   return(0);
}

/******************************************************************************/
void info_field(){

    double totalMass=0.0;
    for(int it=0;it<ntype;it++){
        totalMass+=massit[it]*nionit[it];
    }
    double rho; /* mass density */
    rho = totalMass * simuCell.inveOmega;

    if (ionode){
        SEPARATOR;
        printf("field info\n");
        LSEPARATOR;
        printf("mass of types :\n");
        for(int it=0;it<ntype;it++){
            printf("%d %s                   = %.5f\n",it,atypit[it],massit[it]);
        }
        LSEPARATOR;
        printf("total mass            = %.5f a.m     \n",totalMass);               
        printf("density               = %.5f g/cm^3  \n",rho*rho_unit);
        printf("density(N)            = %.5f ions/A^3\n",rhoN); 
        putchar('\n');
    }

}

/******************************************************************************/
void init_field(char* controlfn){
    
    read_field(controlfn);
    if (lnmlj) {
        lnonbonded=true;
    }
    info_field();
    if (lnmlj) init_nmlj(controlfn);
}

/******************************************************************************/
void engforce()
{
    statime(2);
    if (lnmlj) {
        engforce_nmlj_pbc(&u_lj,&pvir_lj,tau_lj);
    }
    statime(3);
    mestime(&engforceCPUtime,3,2);


    if (lcoul) {
        multipole_ES(qia,dipia,quadia);
    }







}


