#include <stdio.h>
#include <string.h>

#include "constants.h"
#include "config.h"
#include "thermo.h"
#include "field.h"
#include "nmlj.h"
#include "cell.h"
#include "timing.h"
#include "md.h"
#include "io.h"

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
        // mass 
        if (strcmp(buffer,"mass") == 0 ) {
            for(int it=0;it<ntype;it++){
                fscanf(fp,"%lf",&mass[it]);
            }
        } 
        if (strcmp(buffer,"lnmlj") == 0 ) {
            fscanf(fp,"%s",buffer);
            lnmlj=check_FTstring("lnmlj",buffer)?true:false; 
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
    rho = totalMass * simuCell.inveOmega;

    if (ionode){
        SEPARATOR;
        printf("field info\n");
        LSEPARATOR;
        printf("mass of types :\n");
        for(int it=0;it<ntype;it++){
            printf("%d %s                   = %.5f\n",it,atypei[it],mass[it]);
        }
        LSEPARATOR;
        printf("total mass            = %.5f a.m     \n",totalMass);               
        printf("density               = %.5f g/cm^3  \n",rho*rho_unit);
        printf("density(N)            = %.5f ions/A^3\n",rhoN); 
        putchar('\n');
    }

}

void init_field(char* controlfn){
    read_field(controlfn);
    info_field();
    init_nmlj(controlfn);
}


void engforce()
{
    statime(2);
    if (lnmlj) engforce_nmlj_pbc(&u_lj,&vir_lj);
    statime(3);
    mestime(&engforceCPUtime,3,2);
}


