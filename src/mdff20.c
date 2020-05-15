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
    time_t starting_time, finishing_time;
    clock_t t1=clock();
    char* Cstarting_time, *Cfinishing_time;
    
    starting_time = time(NULL);
    Cstarting_time=strdup(asctime(localtime(&starting_time)));

#ifdef MPI
    MPI_Init(NULL,NULL);
    init_io();
    int myrank,numprocs;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
#else
    int myrank=0; int numprocs=1;
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
    headerstdout(Cstarting_time,numprocs);


    //à mettre ailleur
    psimu_cell=&simu_cell;

    //read configuration from file POSFF
    //allocate main quantities when nion is known
    read_posff(psimu_cell);

    // atom decomposition
    patom_dec=&atom_dec;
    do_split(nion,numprocs,myrank,patom_dec,"atoms");
    
//    read_config(controlfn); //needed ?

    // main md parameters
    init_md(controlfn);
    
    init_field(controlfn);
    
    //some initialize quantities à regrouper ailleur
    init_velocities();

    // main md function
    run_md();



    //
    double elapsed_time=((double) clock()-t1)/CLOCKS_PER_SEC;
    printf("\n");
    SEPARATOR;
    finishing_time = time(NULL);
    Cfinishing_time=strdup(asctime(localtime(&finishing_time)));
    
    io_node printf("Elapsed time : %f (s)\n", elapsed_time );
    io_node printf("Date         : %s", Cfinishing_time);

    free_config();
    free(Cstarting_time);
    free(Cfinishing_time);

#ifdef MPI
    MPI_Finalize();
#endif

    return EXIT_SUCCESS;
}
