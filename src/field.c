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
#include "multipole.h"

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
            lcoul=check_boolstring("lcoul",buffer); 
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
                fscanf(fp,"%lf %lf %lf",&dipit[it][0],&dipit[it][1],&dipit[it][2]);
            }
        } 
        // quadrupole on type
        if (strcmp(buffer,"quadit") == 0 ) {
            for(int it=0;it<ntype;it++){
                fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                          &quadit[it][0][0],&quadit[it][0][1],&quadit[it][0][2],
                          &quadit[it][1][0],&quadit[it][1][1],&quadit[it][1][2],
                          &quadit[it][2][0],&quadit[it][2][1],&quadit[it][2][2]);
            }
        } 
        if (strcmp(buffer,"alphaES") == 0 ) {
            fscanf(fp,"%lf",&alphaES);
        } 
        if (strcmp(buffer,"kES") == 0 ) {
            for(int k=0;k<3;k++){
                fscanf(fp,"%d",&kES[k]);
            }
        } 
        if (strcmp(buffer,"epsw") == 0 ) {
            fscanf(fp,"%lf",&epsw);
        } 
        if (strcmp(buffer,"epsw") == 0 ) {
            fscanf(fp,"%lf",&epsw);
        } 
        if (strcmp(buffer,"lautoES") == 0 ) {
            fscanf(fp,"%s",buffer);
            lautoES=check_boolstring("lautoES",buffer);
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
        printf("point charges:\n");
        for (int it=0;it<ntype;it++) {
            printf("q_%s                  = %.5f \n",atypit[it],qit[it]);
            printf("mu_%s                 = ",atypit[it]);
            for(int k=0;k<3;k++){
                printf("%.5f ",atypit[it],dipit[it][k]);
            }
            putchar('\n');
            printf("theta_%s              = \n",atypit[it]);
            for(int j=0;j<3;j++){
                printf("                      ");
                for(int k=0;k<3;k++){
                    printf(" %.5f",quadit[it][j][k]);
                }
                putchar('\n');
            }
            putchar('\n');
        }
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
    if (lautoES) set_autoES();
    if (lcoul) init_multipole();
}

/******************************************************************************/
void engforce()
{
    statime(2);
    if (lnmlj) {
        printf("here NMLJ\n");
        engforce_nmlj_pbc(&u_lj,&pvir_lj,tau_lj);
    }
    statime(3);
    mestime(&engforceCPUtime,3,2);


    if (lcoul) {
        printf("here ES\n");
        multipole_ES(qia,dipia,quadia,&u_coul);
    }







}


