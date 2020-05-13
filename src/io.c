#include<stdio.h>
#include <time.h>
#include "io.h"

void headerstdout(time_t starting_time){

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
    printf("Date :       %s", asctime(localtime(&starting_time)));

    printf("MOLECULAR DYNAMICS ...for fun in 2020 in now in C            \n");
    printf("filipe.manuel.vasconcelos@gmail.com                          \n");

}
