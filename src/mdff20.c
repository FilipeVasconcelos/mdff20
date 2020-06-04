#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef MPI
#include <mpi.h>
#endif
#ifdef OMP 
#include <omp.h>
#endif
#include "color_mdff.h"
#include "constants.h"
#include "global.h"
#include "cell.h"
#include "rand.h"
#include "config.h"
#include "field.h"
#include "nmlj.h"
#include "kinetic.h"
#include "verlet.h"
#include "md.h"
#include "io.h"
#include "timing.h"
#include "thermo.h"
#include "multipole.h"
#include "kspace.h"

/******************************************************************************/
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
#elif OMP
    double startingTime,finishingTime;
    startingTime = omp_get_wtime();
    myrank=0;numprocs=1;
    init_io();
#else
    clock_t startingTime=clock();
    myrank=0;numprocs=1;
    init_io();
#endif
    // init random number generator (velocities)
    //init_rand(startingDate);
    // main constants of the code
    gen_constants();

    if(argc < 2)
    {
        io_node printf(BLU"Usage :"RES" ./mdff20.x <filename>\n");
        exit(-1);
    }
    //main control file
    controlfn=argv[1];

    // header output à la MDFF 
    headerstdout(pstartingDate,numprocs);

    init_global(controlfn);
    //read configuration from file POSFF
    //allocate main quantities when nion is known
    read_config();

    // main md parameters
    init_md(controlfn);
    // main field parameters
    init_field(controlfn);
    init_config();

    // parallelization atom decomposition
    do_split(nion,numprocs,myrank,&atomDec,"atoms");
    do_split(kcoul.nk,numprocs,myrank,&kcoul.kptDec,"k-points");
    // verlet list
    if (lverletL) check_verletlist();

    if (!lstatic) {
        //lets give ions some velocities
        init_velocities();
        // main md function
        run_md();
    }
    else {
        double tempi;
        e_kin=calc_kin();
        io_node printf("static properties\n");
        /* energie, force and pressure */
        engforce();
        if(ionode){
            SEPARATOR;
            printf("properties at t=0\n");
            SEPARATOR;
            putchar('\n');
        }
        info_thermo(0,NULL); 
    }
   
    /* ------------------------------------- */ 
    info_timing();
    finishingDate = time(NULL);
    pfinishingDate=strdup(asctime(localtime(&finishingDate)));
    if (ionode) {
        printf("\n");
        SEPARATOR;
#ifdef MPI
        finishingTime = MPI_Wtime();
        double elapsedTime=finishingTime-startingTime; 
#elif OMP
        finishingTime = omp_get_wtime();
        double elapsedTime=finishingTime-startingTime; 
#else
        double elapsedTime=((double) clock()-startingTime)/CLOCKS_PER_SEC;
#endif
        printf("Elapsed time : %5.3f (s)\n", elapsedTime );
        printf("Date         : %s", pfinishingDate);
    }
    free_config();
    free_md();
    free_multipole();
    free(pstartingDate);
    free(pfinishingDate);
#ifdef MPI
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
