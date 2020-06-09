#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "constants.h"
#include "config.h"
#include "thermo.h"
#include "field.h"
#include "nonbonded.h"
#include "nmlj.h"
#include "bhmftd.h"
#include "ewald.h"
#include "cell.h"
#include "timing.h"
#include "md.h"
#include "io.h"
#include "tools.h"
#include "coulombic.h"

/******************************************************************************/
int read_field(char* controlfn)
{
    int kqch,kdip,kqua,kmass;
    kqch=0;
    kdip=0;
    kqua=0;
    kmass=0;

    char buffer[MAX_LEN+1];
    FILE * fp;
    fp = fopen (controlfn, "r");
    if (NULL == fp )  {
        pError("opening control file");
        return (-1);
    }
    while (EOF != fscanf(fp, "%s\n", buffer)) { 
        if (strcmp(buffer,"lbhmft") == 0 ) {
            fscanf(fp,"%s",buffer);
            lbhmft=check_boolstring("lbhmft",buffer); 
        } 
        if (strcmp(buffer,"lbhmftd") == 0 ) {
            fscanf(fp,"%s",buffer);
            lbhmftd=check_boolstring("lbhmftd",buffer); 
        } 
        if (strcmp(buffer,"lnmlj") == 0 ) {
            fscanf(fp,"%s",buffer);
            lnmlj=check_boolstring("lnmlj",buffer); 
        } 
        if (strcmp(buffer,"lcoulombic") == 0 ) {
            fscanf(fp,"%s",buffer);
            lcoulombic=check_boolstring("lcoulombic",buffer); 
        } 
        // mass of type
        if (strcmp(buffer,"massit") == 0 ) {
            fscanf(fp,"%lf",&massit[kmass]);
            kmass+=1;
        } 
        // charges on type
        if (strcmp(buffer,"qit") == 0 ) {
            fscanf(fp,"%lf",&qit[kqch]);
            kqch+=1;
        } 
        // dipole on type
        if (strcmp(buffer,"dipit") == 0 ) {
            fscanf(fp,"%lf %lf %lf",&dipit[kdip][0],&dipit[kdip][1],&dipit[kdip][2]);
            kdip+=1;
        } 
        // quadrupole on type
        if (strcmp(buffer,"quadit") == 0 ) {
           fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                  &quadit[kqua][0][0],&quadit[kqua][0][1],&quadit[kqua][0][2],
                  &quadit[kqua][1][0],&quadit[kqua][1][1],&quadit[kqua][1][2],
                  &quadit[kqua][2][0],&quadit[kqua][2][1],&quadit[kqua][2][2]);
            kqua+=1;
        } 
        //ewald sum param 
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
    if ( (kmass != ntype) && (kmass !=0) ) {
        pError("massit not found\n");
        printf("kmass %d != ntype %d\n",kmass,ntype);
        exit(-1);
    } 
    if ( (kqch != ntype) && (kqch !=0) ) {
        pError("qit not found\n");
        printf("kdip %d != ntype %d\n",kqch,ntype);
        exit(-1);
    } 
    if ( (kdip != ntype) && (kdip !=0) ) {
        pError("dipit not found\n");
        printf("kdip %d != ntype %d\n",kdip,ntype);
        exit(-1);
    } 
    if ( (kqua != ntype) && (kqua !=0) ) {
        pError("quadit not found\n");
        printf("kdip %d != ntype %d\n",kqua,ntype);
        exit(-1);
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

    lqch = false; ldip= false; lqua = false;
    for(int it=0; it< ntype ; it++){
        if (qit[it] != 0.0 ) lqch = true;
        for (int j=0;j<3;j++){
            if (dipit[it][j] != 0.0 ) ldip = true;
            for(int k=0;k<3;k++){
                if (quadit[it][j][k] != 0.0 ) lqua = true;
            }
        }
    } 
    if ( ( lcoulombic ) && ! ( ( lqch)  || ( ldip ) || ( lqua ) ) ) {
        pError("with lcoulombic :  charges, dipoles or quadrupoles need to be set\n");
        exit(-1);
    }


    if (ionode){
        SEPARATOR;
        printf("field info\n");
        LSEPARATOR;
        printf("mass of types :\n");
        for(int it=0;it<ntype;it++){
            printf("m_%-2s                  = %.5f a.m\n",atypit[it],massit[it]);
        }
        LSEPARATOR;
        printf("total mass            = %.5f a.m     \n",totalMass);               
        printf("density               = %.5f g/cm^3  \n",rho*rho_unit);
        printf("density(N)            = %.5f ions/A^3\n",rhoN); 
        putchar('\n');
        if (lqch) {
            printf("point charges :\n");
            printf("--------------\n");
            for (int it=0;it<ntype;it++) {
                printf("q_%-2s                  = %8.5f \n",atypit[it],qit[it]);
            }
            putchar('\n');
        }
        if (ldip) {
            printf("dipoles :\n");
            printf("--------\n");
            for (int it=0;it<ntype;it++) {
                printf("mu_%s                 = ",atypit[it]);
                for(int k=0;k<3;k++){
                    printf("%8.5f ",dipit[it][k]);
                }
                putchar('\n');
            }
            putchar('\n');
        }
        if (lqua) {
            printf("quadrupoles :\n");
            printf("------------\n");
            for (int it=0;it<ntype;it++) {
                printf("theta_%s              = \n",atypit[it]);
                for(int j=0;j<3;j++){
                    printf("                      ");
                    for(int k=0;k<3;k++){
                        printf(" %8.5f",quadit[it][j][k]);
                    }
                    putchar('\n');
                }
                putchar('\n');
            }
        }
    }
}


/******************************************************************************/
void init_field(char* controlfn){
    
    read_field(controlfn);

    info_field();

    /* NON BONDED POTENTIAL */
    if ( (lnmlj) || (lbhmft) || (lbhmftd) ) {
        srcutsq=cutshortrange*cutshortrange;
        lnonbonded=true;
        init_nonbonded(controlfn);
        if (lnmlj) init_nmlj(controlfn);
        if ( (lbhmft) || (lbhmftd) ) init_bhmftd(controlfn);
    }
    /* ELECTROSTATIC POTENTIAL */
    if (lcoulombic) {
        lrcutsq=cutlongrange*cutlongrange;
        if (lautoES) set_autoES();
        init_coulombic();
    }
}

/******************************************************************************/
void engforce()
{
    for(int ia=0;ia<nion;ia++){
        fx[ia]=0.0;fy[ia]=0.0;fz[ia]=0.0;
    }
    statime(2);
    if (lnmlj) {
        engforce_nmlj_pbc(&u_lj,&pvir_lj,tau_lj);
    }
    statime(15);
    mestime(&engforce_nmljCPUtime,15,2);
    if ( (lbhmft) || (lbhmftd) ) {
        engforce_bhmftd_pbc(&u_bhmftd,&pvir_bhmftd,tau_bhmftd);
    }
    statime(15);
    mestime(&engforce_bhmftdCPUtime,15,2);

    if (lcoulombic) {
        double (*ef)[3];
        double (*efg)[3][3];
        ef=malloc(nion*sizeof(*ef));
        efg=malloc(nion*sizeof(*efg));
        multipole_ES(qia,dipia,quadia,&u_coul,&pvir_coul,tau_coul,ef,efg);
        free(ef);
        free(efg);
    }
    statime(3);
    mestime(&engforceCPUtime,3,2);
    mestime(&engforce_coulCPUtime,3,15);

}


