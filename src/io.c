#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/utsname.h>
#include <time.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "io.h"

/******************************************************************************/
void headerstdout(char *starting_time, int numprocs){

    char *username = getenv("USER");
    char *plural;
    plural = ((numprocs>1) ? "s" : "\0");
    struct utsname hostname;
    uname(&hostname);

    if (ionode) {
      printf("                            \\\\|//                          \n");
      printf("                           -(o o)-                           \n");
      printf("========================oOO==(_)==OOo========================\n");
      printf("     ____    ____ ______  ________ ________                  \n");
      printf("    |_   \\  /   _|_   _ `|_   __  |_   __  |                \n");
      printf("      |   \\/   |   | | `. \\| |_ \\_| | |_ \\_|____    ____ \n");
      printf("      | |\\  /| |   | |  | ||  _|    |  _| / ___ `..'    '.  \n");
      printf("     _| |_\\/_| |_ _| |_.' _| |_    _| |_ |_/___) |  .--.  | \n");
      printf("    |_____||_____|______.|_____|  |_____| .'____.| |    | |  \n");
      printf("                                         / /_____|  `--'  |  \n");
      printf("                                         |_______|'.____.'   \n");
      printf("\n");
      SEPARATOR;
      printf("MOLECULAR DYNAMICS ...for fun in 2020 ... but now in C ;)    \n");
      printf("filipe.manuel.vasconcelos@gmail.com                          \n");
      printf("Running on : %d node%s\n",numprocs,plural);
      printf("by user    : %s\n",username);
      printf("host       : %s %s %s\n",hostname.nodename,
                                       hostname.sysname,
                                       hostname.machine);
      printf("Date       : %s", starting_time);
    }
}
/******************************************************************************/
void init_io(){
    int myrank =0;
#ifdef MPI
    int numprocs;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
#endif
    if (myrank==0) {
        ionode = true;
    }
    else{
        ionode = false;
    }
}

/******************************************************************************/
int iopnode(int step,int npas,int nprint){
//    printf("step : %d %d %d %d\n",step,npas,nprint,step%nprint==0);
    return ( (step % nprint==0 || step == npas ) && (step > 0) );
}


