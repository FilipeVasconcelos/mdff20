#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

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
    time_t starting_time, finishing_time;
    clock_t t1=clock();
    char* Cstarting_time, *Cfinishing_time;

#ifdef MPI_
    MPI_Init(NULL,NULL);

    int myrank,numprocs;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    printf("proc : %d numprocs : %d\n",myrank,numprocs);
#else
    int numprocs=1;
#endif

    init_rand();
    if(argc < 2)
    {
        printf("Usage : ./mdff20 <filename>\n");
        exit(0);
    }
    controlfn=argv[1];
    starting_time = time(NULL);
    Cstarting_time=strdup(asctime(localtime(&starting_time)));
   
    // header output à la MDFF 
    headerstdout(Cstarting_time,numprocs);

    // main constants of the code
    gen_constants();

    //à mettre ailleur
    psimu_cell=&simu_cell;

    //read configuration from file POSFF
    //allocate main quantities when nion is known
    read_posff(psimu_cell);

    //à mettre dans POSFF


    
//    read_config(controlfn); //needed ?

    //read md parameters 
    read_md(controlfn);
    //read field parameters 
    read_field(controlfn);

    //some initialize quantities à regrouper ailleur
    init_nmlj();
    init_velocities();
    maxwellboltzmann_velocities();

    // main md function
    run_md();
    

    //
    double elapsed_time=((double) clock()-t1)/CLOCKS_PER_SEC;
    printf("\n");
    SEPARATOR;
    finishing_time = time(NULL);
    Cfinishing_time=strdup(asctime(localtime(&finishing_time)));
    printf("Elapsed time : %f (s)\n", elapsed_time );
    printf("Date         : %s", Cfinishing_time);

    free(vx);
    free(vy);
    free(vz);
    free(massia);
    //free(atype); double free corruption ???
    free(Cstarting_time);
    free(Cfinishing_time);

#ifdef MPI_
    MPI_Finalize();
#endif

    return EXIT_SUCCESS;
}
