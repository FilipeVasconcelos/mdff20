#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "color_mdff.h"
#include "constants.h"
#include "global.h"
#include "io.h"
#include "verlet.h"

int read_global(char* controlfn)
{
   char buffer[MAX_LEN+1];
   FILE * fp;
   fp = fopen (controlfn, "r");
   if (NULL == fp )  {
       perror("opening database file");
       return (-1);
   }
   while (EOF != fscanf(fp, "%s\n", buffer)) { 
        if (strcmp(buffer,"lverletL") == 0 ) {
            fscanf(fp,"%s",buffer);
            lverletL=check_FTstring("lverletL",buffer)?true:false; 
        } 
        if (strcmp(buffer,"cutshortrange") == 0 ) {
            fscanf(fp,"%lf",&cutshortrange);
        } 
        if (strcmp(buffer,"skindiff") == 0 ) {
            fscanf(fp,"%lf",&skindiff);
        } 
   }
   fclose(fp);
   return(0);
}

void info_global(){
    if (ionode){
        SEPARATOR;
        printf("global info\n");
        LSEPARATOR;
        printf("verlet list           = %s \n",lverletL?"true":"false");               
        printf("cutshortrange         = %-6.2f \n",cutshortrange);
        printf("skindiff              = %-6.2f \n",skindiff);
        putchar('\n');
    }
}

void init_global(char* controlfn){
    read_global(controlfn);
    info_global();
    if (lverletL) { 
        verlet_nb=allocate_verletlist("vnlnb");
        verlet_nb->cut=cutshortrange; 
        printf("test %d\n",verlet_nb->point[0]);
    }
}


