#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "constants.h"
#include "cell.h"
#include "rand.h"
#include "read_posff.h"
#include "config.h"
#include "field.h"
#include "nmlj.h"
#include "kinetic.h"
#include "md.h"
#include "io.h"


int main(int argc, char *argv[])
{
    time_t startingDate, finishingDate;
    char *pstartingDate, *pfinishingDate;
    startingDate = time(NULL);
    pstartingDate=strdup(asctime(localtime(&startingDate)));

#ifdef MPI
    MPI_Init(NULL,NULL);
    double startingTime, finishingTime; 
    startingTime = MPI_Wtime(); 
    init_io();
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
#else
    clock_t startingTime=clock();
    myrank=0;numprocs=1;
    init_io();
#endif
    // init random number generator (velocities)
    init_rand();
    // main constants of the code
    gen_constants();

    if(argc < 2)
    {
        io_node printf("Usage : ./mdff20 <filename>\n");
        exit(0);
    }
    controlfn=argv[1];

    // header output à la MDFF 
    headerstdout(pstartingDate,numprocs);

    //à mettre ailleur
    psimuCell=&simuCell;

    //read configuration from file POSFF
    //allocate main quantities when nion is known
    read_posff(psimuCell);

    // atom decomposition
    patomDec=&atomDec;
    do_split(nion,numprocs,myrank,patomDec,"atoms");
    
//    read_config(controlfn); //needed ?

    // main md parameters
    init_md(controlfn);
    
    init_field(controlfn);
    
    //some initialize quantities à regrouper ailleur
    //init_velocities();

    // main md function
    run_md();
    
    
    finishingDate = time(NULL);
    pfinishingDate=strdup(asctime(localtime(&finishingDate)));

    if (ionode) {
        printf("\n");
        SEPARATOR;
#ifdef MPI
        finishingTime = MPI_Wtime();
        double elapsedTime=finishingTime-startingTime; 
#else
        double elapsedTime=((double) clock()-startingTime)/CLOCKS_PER_SEC;
#endif
        printf("Elapsed time : %5.3f (s)\n", elapsedTime );
        printf("Date         : %s", pfinishingDate);
    }
    free_config();
    free(pstartingDate);
    free(pfinishingDate);
#ifdef MPI
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
