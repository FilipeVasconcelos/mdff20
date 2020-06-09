#include <stdio.h>
#include <string.h>
#include "io.h"
#include "tools.h"
#include "pim.h"

/******************************************************************************/
int read_pim (char * controlfn) 
{
    char buffer[MAX_LEN+1];
    FILE * fp;

    fp = fopen (controlfn, "r");

    if (NULL == fp ) {
        pError("opening control file (reading pim) ");
        return (-1);
    }

    while (EOF != fscanf(fp, "%s\n", buffer)) {
        if (strcmp(buffer,"algo_pim") == 0 ) {
            fscanf(fp,"%s",buffer);
            algo_pim=check_string("algo_pim",buffer,allwd_algo_pim, ALLWD_ALGO_PIM_STR); 
        } 
        if (strcmp(buffer,"algo_extrapolate_dipole") == 0 ) {
            fscanf(fp,"%s",buffer);
            algo_extrapolate_dipole=check_string("algo_extrapolate_dipole",buffer,allwd_algo_extrapolate_dipole,
                                                           ALLWD_ALGO_EXTRAPOLATE_DIPOLE_STR); 
        } 
        if (strcmp(buffer,"conv_tol_ind") == 0 ) {
            fscanf(fp,"%lf",&conv_tol_ind);
        } 
        if (strcmp(buffer,"extrapolate_order") == 0 ) {
            fscanf(fp,"%d",&extrapolate_order);
        } 
   }
   fclose(fp);
   return 0;
}


void check_pim(){

}

/******************************************************************************/
void init_pim(char * controlfn){
    /* gen allowed input strings for algoPIM */
    strcpy(allwd_algo_pim[0],"scf");
    strcpy(allwd_algo_pim[1],"scfKO");
    /* gen allowed input strings for algo_ext_dipole */
    strcpy(allwd_algo_extrapolate_dipole[0],"poly");
    strcpy(allwd_algo_extrapolate_dipole[1],"aspc");

    info_pim();
}

/******************************************************************************/
void info_pim(){
    if (ionode) {
        SEPARATOR;
        printf("pim info \n");
        LSEPARATOR;
        putchar('\n');
    }
}


void momentpolaSCF(double (*mu)[3]){

}
