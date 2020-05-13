#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"
#include "config.h"

int read_config (char* controlfn) {
   char buffer[10];
   int data;
   FILE * fp;


   fp = fopen (controlfn, "r");
   
   if (NULL == fp ) 
   {
       perror("opening input file");
       return (-1);
   }
   /*
   while (EOF != fscanf(fp, "%s %d\n", buffer,&data))
   {
       if (strcmp(buffer,"nion") == 0 )
       {
           nion=data;
           printf("number of ions  : %d\n", nion );
       } 

        if (strcmp(buffer,"type") == 0 )
        {
            ntype=data;
            printf("number of types : %d\n",ntype );
        } 
   }
    */

   fclose(fp);
   
   return 0;
}


void init_config()
{

    int ccs,cc;

    itype=malloc(nion*sizeof*itype);


    printf("inside init_config %d %d\n",nion,ntype);
    cc=0;
    for ( int it=0; it <= ntype-1;it++)
    {
        ccs = cc;
        cc  = cc + natmi[it];
      //  printf("it : %d %d %d %.2f\n",it,ccs, cc,mass[it]);

        for ( int ia=ccs;ia<=(cc-1);ia++)
        {
            massia[ia] = mass[it];
            invemassia[ia] = 1.0/mass[it];
            itype[ia]  = it;
    //        printf("ia : %d %d %d %d %.2f %.2f \n",it,ia,ccs,cc,mass[it],massia[ia]);
//            printf("it : %d %d %d\n",it,ia,ccs, cc,mass[it],massia[ia]);
        }

    }

}
