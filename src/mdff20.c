#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "constants.h"
#include "cell.h"
#include "rand.h"
#include "read_posff.h"
#include "config.h"
#include "field.h"
#include "nmlj.h"
#include "kinetic.h"
#include "md.h"


int main(int argc, char *argv[])
{
    time_t starting_time, finishing_time;
    clock_t t1=clock();
    double elapsed;
    char * controlfn;

    if(argc < 2)
    {
        printf("Usage : ./mdff20 <filename>\n");
        exit(0);
    }
    controlfn=argv[1];

    starting_time = time(NULL);
    printf("%s\n","===========================================");
    printf("%s\n","                MDFF C2020                 ");
    printf("%s\n","===========================================");
    printf("Date :       %s", asctime(localtime(&starting_time)));
    printf("Input file : %s\n",controlfn);

    gen_constants();

    psimu_cell=&simu_cell;
    read_posff(psimu_cell);
    lattice(psimu_cell);
//    read_config(controlfn);
    read_md(controlfn);
    read_field(controlfn);

    init_nmlj();
    init_config();
    /*
    printf("%f\n",simu_cell.B[0][0]);
    printf("%f\n",simu_cell.A[0][0]);
    exit(1);
    */
    init_velocities();
    //print_velocities();
    maxwellboltzmann_velocities();
    //print_velocities();

    run_md();
    
    finishing_time = time(NULL);
    elapsed = starting_time - finishing_time;
    printf("Elapsed time : %f\n", ((double) clock()-t1)/CLOCKS_PER_SEC );
    printf("Date :       %s", asctime(localtime(&finishing_time)));

    free(vx);
    free(vy);
    free(vz);
    free(massia);
    //free(atype);

    return EXIT_SUCCESS;
}
