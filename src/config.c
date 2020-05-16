#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "color_mdff.h"
#include "config.h"
#include "constants.h"
#include "cell.h"
#include "tools.h"
#include "io.h"

//-----------------------------------------------------------------------------
int read_config (char* controlfn) {

//   char buffer[10];
//   int data;
   FILE * fp;

   fp = fopen (controlfn, "r");
   if (NULL == fp ) {
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


//-----------------------------------------------------------------------------
void alloc_config(){

    onenion=1.0/( (double) nion) ;
    onethirdnion=1.0/( 3.0*(double) nion) ;
    massia= malloc(nion*sizeof *massia );
    invemassia= malloc(nion*sizeof *invemassia );
    atype= malloc(nion*sizeof *atype);
    if (atype==NULL){
        io_node printf("**atype malloc failure\n");
    }
    for (int ia=0;ia<nion;ia++)  {
        atype[ia]=NULL;
        atype[ia]=malloc((MAX_LEN+1)*sizeof(*atype));
        if (atype[ia] == NULL){
            io_node printf("atype %d malloc failure\n",ia);
        }
    }
    vx= malloc(nion*sizeof(*vx));
    vy= malloc(nion*sizeof(*vy));
    vz= malloc(nion*sizeof(*vz));
    rx= malloc(nion*sizeof(*rx));
    ry= malloc(nion*sizeof(*ry));
    rz= malloc(nion*sizeof(*rz));
    fx= malloc(nion*sizeof(*fx));
    fy= malloc(nion*sizeof(*fy));
    fz= malloc(nion*sizeof(*fz));
    rxs= malloc(nion*sizeof(*rxs));
    rys= malloc(nion*sizeof(*rys));
    rzs= malloc(nion*sizeof(*rzs));
    for(int ia=0; ia<nion; ia++) {
        massia[ia]=1.0;
        invemassia[ia]=1.0;
        rx[ia]=0.0;ry[ia]=0.0;rz[ia]=0.0;
        vx[ia]=0.0;vy[ia]=0.0;vz[ia]=0.0;
        fx[ia]=0.0;fy[ia]=0.0;fz[ia]=0.0;
        rxs[ia]=0.0;rys[ia]=0.0;rzs[ia]=0.0;
    }
    init_config();
}

//-----------------------------------------------------------------------------
void free_config(){
    free(massia);
    free(invemassia);
    for (int ia=0;ia<nion;ia++)  {
        free(atype[ia]);
    }
    free(atype);
    free(rx);free(ry);free(rz);
    free(vx);free(vy);free(vz);
    free(fx);free(fy);free(fz);
    free(rxs);free(rys);free(rzs);
}

bool config_input_ok(){
    
    //check if sum(natmi) == nion
//    int sumnatmi=0;
//    for(int it=0;it<ntype;it++){
//        sumnatmi+=natmi[it];
//    }
    int sumnatmi=sum(natmi,ntype);
    if ( sumnatmi != nion) {
        pError("natmi and nions does not correspond\n");
        printf("%d %d\n",sumnatmi,nion);
        return false;
    }
    return true;
}

//-----------------------------------------------------------------------------
void init_config()
{
    if ( ! config_input_ok() ) {
        exit(1);
    }
    int ccs,cc;
    itype=malloc(nion*sizeof*itype);
    cc=0;
    for ( int it=0; it <= ntype-1;it++) {
        ccs = cc;
        cc  = cc + natmi[it];
        for ( int ia=ccs;ia<=(cc-1);ia++) {
            massia[ia] = mass[it];
            invemassia[ia] = 1.0/mass[it];
            itype[ia]  = it;
        }
    }
    info_config();

}
//-----------------------------------------------------------------------------
void info_config(){

    if (ionode) {  
        printf("configname                      : %s\n",configname);
        printf("number of types                 : %d\n",ntype);
        printf("number of ions                  : %d\n",nion);
        putchar('\n');
    }
    info_simuCell();
}

