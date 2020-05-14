#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "constants.h"
#include "config.h"
#include "cell.h"

int read_posff(CELL *Cell)
{

    char buffer[MAX_LEN+1];
    char cpos[MAX_LEN+1];
    int data;
    double f1,f2,f3,f4,f5,f6,f7,f8,f9;
    FILE * fp;

    fp = fopen ("POSFF","r");

    if (NULL == fp ) {
        perror("opening database file\n");
        return (-1);
    }
    // print out info to stdout
    printf("reading configuration\n");
    printf("config from file POSFF\n");
    
    // reading nion number of ions in POSFF 
    fscanf(fp, "%d", &nion);

    
    // reading config name
    fscanf(fp, "%s", buffer);
    strcpy(configname,buffer); 

    // reading cell parameters
    fscanf(fp, "%lf %lf %lf\n",&f1,&f2,&f3);
    fscanf(fp, "%lf %lf %lf\n",&f4,&f5,&f6);
    fscanf(fp, "%lf %lf %lf\n",&f7,&f8,&f9);
    Cell->A[0][0]=f1;   Cell->A[0][1]=f2;   Cell->A[0][2]=f3;
    Cell->A[1][0]=f4;   Cell->A[1][1]=f5;   Cell->A[1][2]=f6;
    Cell->A[2][0]=f7;   Cell->A[2][1]=f8;   Cell->A[2][2]=f9;
    lattice(Cell);

    // ------------------
    // type informations 
    // ------------------
    //
    //ntype : number of types in POSFF
    fscanf(fp, "%d", &ntype); //ntype

    //can be call only when nion and ntype in known
    alloc_config();

    printf("here !!!!\n");
    //atypei : ion types char 
    for (int it=0;it<ntype;it++) {
        fscanf(fp,"%s",buffer);
        strcpy(atypei[it],buffer);
        printf("%d %s\n",it,atypei[it]);
    }
    //natmi : ions per type
    for (int it=0;it<ntype;it++) {
        fscanf(fp,"%d",&data);
        natmi[it]=data;
        printf("%d %d\n",it,natmi[it]);
    }

    // Direct or Cartesian
    fscanf(fp, "%s", buffer);
    strcpy(cpos,buffer);

    // reading positions of ions 
    for (int ia=0;ia<nion;ia++) {
        fscanf(fp,"%s %lf %lf %lf",buffer,&f1,&f2,&f3);
        strcpy(atype[ia],buffer);
        rx[ia]=f1;ry[ia]=f2;rz[ia]=f3;
    }

    if (strcmp(cpos,"Direct") == 0 ) {
        printf("Direct to Cartesian\n");
        dirkar(nion, rx, ry, rz , Cell->A);
    }

    //closing POSFF
    if (fclose(fp))     { 
       printf("error closing file."); 
       exit(-1); 
    }
   /* 
    printf("configname  : %s\n", buffer);
    printf("number of types  : %d\n", ntype );

    printf("number of ions  : %d\n", nion );
    printf("leaving read_posfff\n");
    */
    return 0;
}
