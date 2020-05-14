#include <stdio.h>
#include <stdlib.h>
#include <sys/utsname.h>
#include <time.h>
#include "io.h"

void headerstdout(time_t starting_time, int numprocs){

    char *username = getenv("USER");
    char *plural;
    plural = ((numprocs>1) ? "s" : "\0");
    struct utsname hostname;
    uname(&hostname);

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
    printf("Runnin on : %d node%s\n",numprocs,plural);
    printf("by user   : %s\n",username);
    printf("host      : %s %s %s\n",hostname.nodename,hostname.sysname,hostname.machine);
    printf("Date :       %s", asctime(localtime(&starting_time)));

}
