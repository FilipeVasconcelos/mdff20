#include<stdio.h>
#include <string.h>

#include "constants.h"
#include "config.h"
#include "thermo.h"
#include "field.h"
#include "nmlj.h"
#include "cell.h"
#include "io.h"

int read_field(char* controlfn)
{
   char buffer[MAX_LEN+1];
   double data;
   FILE * fp;
   fp = fopen (controlfn, "r");
   if (NULL == fp )  {
       perror("opening database file");
       return (-1);
   }
   while ( (EOF != fscanf(fp, "%s %lf\n", buffer,&data)) ||
           (EOF != fscanf(fp, "%s\n", buffer)) )  {
        if (strcmp(buffer,"cutshortrange") == 0 ) {
            cutshortrange=data;
        } 
        if (strcmp(buffer,"trunctype") == 0 ) {
            trunctype=data;
        } 
        if (strcmp(buffer,"skindiff") == 0 ) {
            skindiff=data;
        } 
        if (strcmp(buffer,"mass") == 0 ) {
            mass[0]=data;
        } 
        if (strcmp(buffer,"lnmlj") == 0 ) {
            lnmlj=true;
        } 
   }
   fclose(fp);
   return(0);
}

void info_field(){

    double totalMass=0.0;
    for(int it=0;it<ntype;it++){
        totalMass+=mass[it]*natmi[it];
    }
    double rho; /* mass density */
    rho = totalMass * simu_cell.inveOmega;

    if (ionode){
        printf("total mass            = %12.4f a.m     \n",totalMass);               
        printf("density               = %12.4f g/cm^3  \n",rho*rho_unit);
        printf("density(N)            = %12.4f ions/A^3\n",rhoN); 
    }

}

void init_field(char* controlfn){
    read_field(controlfn);
    init_nmlj();
}


void engforce()
{
    double ulj;
    engforce_nmlj_pbc(&ulj);
    u_lj=ulj;
}


