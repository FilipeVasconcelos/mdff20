#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef MPI
#include <mpi.h>
#endif

#include "constants.h"
#include "config.h"
#include "cell.h"
#include "io.h"

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
#ifdef MPI
        MPI_Finalize();
#endif
        return (-1);
    }
    // print out info to stdout
    io_node printf("reading configuration from file POSFF\n");
    
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
    rhoN=(double)nion * Cell->inveOmega;

    // ------------------
    // type informations 
    // ------------------
    //
    //ntype : number of types in POSFF
    fscanf(fp, "%d", &ntype); //ntype

    //atypei : ion types char 
    io_node printf("found type information on POSFF :");
    for (int it=0;it<ntype;it++) {
        fscanf(fp,"%s",buffer);
        strcpy(atypei[it],buffer);
        io_node printf("%5s",atypei[it]);
    }
    //natmi : ions per type
    io_node printf("\n                                 ");
    for (int it=0;it<ntype;it++) {
        fscanf(fp,"%d",&data);
        natmi[it]=data;
        io_node printf("%5d",natmi[it]);
    }
    io_node putchar('\n');
    //can be call only when nion and ntype 
    //are known and  natmi,atypei readed
    alloc_config();

    // Direct or Cartesian
    fscanf(fp, "%s", buffer);
    strcpy(cpos,buffer);
    // reading positions of ions 
    for (int ia=0;ia<nion;ia++) {
        fscanf(fp,"%s %lf %lf %lf",buffer,&f1,&f2,&f3);
        strcpy(atype[ia],buffer);
        rx[ia]=f1;ry[ia]=f2;rz[ia]=f3;
    }
    if ((strcmp(cpos,"Direct") == 0 ) || (strcmp(cpos,"D") == 0 )) {
        io_node printf("\natomic positions in Direct coordinates\n");
        dirkar(nion, rx, ry, rz , Cell->A);
    }
    else{
        io_node printf("\natomic positions in Cartesian coordinates\n");
    }
    io_node putchar('\n');

    //closing POSFF
    if (fclose(fp))     { 
       io_node printf("error closing file."); 
#ifdef MPI
        MPI_Finalize();
#endif
       return -1; 
    }
    return 0;
}
