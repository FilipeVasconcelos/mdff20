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

int read_posff()
{

    char cpos[MAX_LEN+1];
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
    fscanf(fp, "%s", configname);
    // reading cell parameters
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            fscanf(fp, "%lf",&simuCell.A[i][j]);
        }
    }
    lattice(&simuCell);
    rhoN=(double)nion * simuCell.inveOmega;

    // ------------------
    // type informations 
    // ------------------
    //ntype : number of types in POSFF
    fscanf(fp, "%d", &ntype); //ntype

    //atypei : ion types char 
    io_node printf("found type information on POSFF :");
    for (int it=0;it<ntype;it++) {
        fscanf(fp,"%s",atypei[it]);
        io_node printf("%5s",atypei[it]);
    }
    //natmi : ions per type
    io_node printf("\n                                 ");
    for (int it=0;it<ntype;it++) {
        fscanf(fp,"%d",&natmi[it]);
        io_node printf("%5d",natmi[it]);
    }
    io_node putchar('\n');
    //can be call only when nion and ntype 
    //are known and  natmi,atypei readed
    alloc_config();

    // Direct or Cartesian
    fscanf(fp, "%s", cpos);
    //strcpy(cpos,buffer);

    // reading positions of ions 
    for (int ia=0;ia<nion;ia++) {
        fscanf(fp,"%s %lf %lf %lf",atype[ia],&rx[ia],&ry[ia],&rz[ia]);
    }
    // if position are in direct coordinates => cartesian coordinates
    if ((strcmp(cpos,"Direct") == 0 ) || (strcmp(cpos,"D") == 0 )) {
        io_node printf("\natomic positions in Direct coordinates\n");
        dirkar(nion, rx, ry, rz , simuCell.A);
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
