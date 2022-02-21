#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "color_mdff.h"
#include "constants.h"
#include "global.h"
#include "io.h"
#include "verlet.h"
#include "tools.h"

/******************************************************************************/
int read_global(char* controlfn)
{
    int c;
    char buffer[MAX_LEN+1];
    FILE * fp;
    fp = fopen (controlfn, "r");
    if (NULL == fp ) {
        pError("opening control file");
        return (-1);
    }
    while (EOF != fscanf(fp, "%s\n", buffer)) {
        if (strcmp(buffer,"lverletL") == 0 ) {
            c=fscanf(fp,"%s",buffer);
            lverletL=check_boolstring("lverletL",buffer);
        }
        if (strcmp(buffer,"lstatic") == 0 ) {
            c=fscanf(fp,"%s",buffer);
            lstatic=check_boolstring("lstatic",buffer);
        }
        if (strcmp(buffer,"lrdf") == 0 ) {
            c=fscanf(fp,"%s",buffer);
            lrdf=check_boolstring("lrdf",buffer);
        }
        if (strcmp(buffer,"lpstress") == 0 ) {
            c=fscanf(fp,"%s",buffer);
            lpstress=check_boolstring("lpstress",buffer);
        }
        if (strcmp(buffer,"cutshortrange") == 0 ) {
            c=fscanf(fp,"%lf",&cutshortrange);
        }
        if (strcmp(buffer,"cutlongrange") == 0 ) {
            c=fscanf(fp,"%lf",&cutlongrange);
        }
        if (strcmp(buffer,"Fposff") == 0 ) {
            c=fscanf(fp,"%s",buffer);
            Fposff=check_string("Fposff",buffer,allwd_Fposff,ALLWD_FORMAT_POSFF_STR);
        }
        if (strcmp(buffer,"skindiff") == 0 ) {
            c=fscanf(fp,"%lf",&skindiff);
        }
        if (strcmp(buffer,"lreduced") == 0 ) {
            c=fscanf(fp,"%s",buffer);
            lreduced=check_boolstring("lreduced",buffer);
        }
    }
    fclose(fp);
    return(0);
}

/******************************************************************************/
void info_global(){

    if (ionode){
        SEPARATOR;
        printf("global info\n");
        LSEPARATOR;
        printf("static calculation    = %s \n",lstatic?"true":"false");
        printf("radial distribution   = %s \n",lrdf?"true":"false");
        printf("stress tensor         = %s \n",lpstress?"true":"false");
        printf("verlet list           = %s \n",lverletL?"true":"false");
        printf("reduced units         = %s \n",lreduced?"true":"false");
        printf("cutshortrange         = %-6.2f \n",cutshortrange);
        printf("cutlongrange          = %-6.2f \n",cutlongrange);
        printf("skindiff              = %-6.2f \n",skindiff);
        printf("Fposff                = %d     \n",Fposff);
        putchar('\n');
    }
}

/******************************************************************************/
void default_global(){
    Fposff=2;
    /* gen allowed input strings for integrator */
    strcpy(allwd_Fposff[0],"rnn"); /* only positions */
    strcpy(allwd_Fposff[1],"rvn"); /* positions + velocities */
    strcpy(allwd_Fposff[2],"rvf"); /* positions + velocities + forces */
    skindiff=0.15;
    lstatic=false;
    lpstress=false;
    lrdf=false;
    lverletL=true;
}

/******************************************************************************/
void check_global(){
    if (lreduced) reduced_units();
}

/******************************************************************************/
void init_global(char* controlfn){
    /* default values */
    default_global();
    read_global(controlfn);
    check_global();
    info_global();
}


