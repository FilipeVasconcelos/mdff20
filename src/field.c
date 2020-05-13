#include<stdio.h>
#include <string.h>

#include "constants.h"
#include "thermo.h"
#include "field.h"
#include "nmlj.h"

int read_field(char* controlfn)
{
   char buffer[MAX_LEN+1];
   double data;
   FILE * fp;

   fp = fopen (controlfn, "r");

   if (NULL == fp ) 
   {
       perror("opening database file");
       return (-1);
   }

   while (EOF != fscanf(fp, "%s %lf\n", buffer,&data))
   {
        if (strcmp(buffer,"cutshortrange") == 0 )
        {
            cutshortrange=data;
            printf("cutshortrange   : %.2f\n", cutshortrange );
        } 
        if (strcmp(buffer,"skindiff") == 0 )
        {
            skindiff=data;
            printf("skindiff        : %.2f\n", skindiff );
        } 
        if (strcmp(buffer,"mass") == 0 )
        {
            mass[0]=data;
            printf("skindiff        : %.2f\n", skindiff );
        } 
   }

   fclose(fp);
   
   return(0);
}


void engforce()
{
    double ulj;

    engforce_nmlj_pbc(&ulj);

    u_lj=ulj;
}
