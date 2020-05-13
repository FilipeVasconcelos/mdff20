#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"
#include "config.h"
#include "cell.h"

int read_posff(CELL *Cell)
{
    char buffer[MAX_LEN+1];
    char cpos[MAX_LEN+1];
    int data;
    double f1,f2,f3,f4,f5,f6,f7,f8,f9;
    FILE * fp;
//   CELL Cell;

    fp = fopen ("POSFF","r");

    if (NULL == fp ) {
        perror("opening database file");
        return (-1);
    }
    fscanf(fp, "%d", &nion);
    printf("number of ions  : %d\n", nion );

    onenion=1.0/( (double) nion) ;
    massia= malloc(nion*sizeof *massia );
    invemassia= malloc(nion*sizeof *invemassia );
    atype= malloc(nion*sizeof **atype);
    for (int ia=0;ia<nion;ia++)  {
        atype[ia]=NULL;
        atype[ia]=malloc((MAX_LEN+1)*sizeof(char*));
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

    fscanf(fp, "%s", buffer);
    printf("system  : %s\n", buffer);
    fscanf(fp, "%lf %lf %lf\n",&f1,&f2,&f3);
    fscanf(fp, "%lf %lf %lf\n",&f4,&f5,&f6);
    fscanf(fp, "%lf %lf %lf\n",&f7,&f8,&f9);
   
    printf("%f %f %f\n",f1,f2,f3);
    printf("%f %f %f\n",f4,f5,f6);
    printf("%f %f %f\n",f7,f8,f9);
   
    Cell->A[0][0]=f1;   Cell->A[0][1]=f2;   Cell->A[0][2]=f3;
    Cell->A[1][0]=f4;   Cell->A[1][1]=f5;   Cell->A[1][2]=f6;
    Cell->A[2][0]=f7;   Cell->A[2][1]=f8;   Cell->A[2][2]=f9;
//   simu_cell=Cell;
 
 /*  
   printf("%lf %lf %lf\n",simu_cell.A[0][0],simu_cell.A[1][0],simu_cell.A[2][0]);
   printf("%lf %lf %lf\n",simu_cell.A[0][1],simu_cell.A[1][1],simu_cell.A[2][1]);
   printf("%lf %lf %lf\n",simu_cell.A[0][2],simu_cell.A[1][2],simu_cell.A[2][2]);
*/
    fscanf(fp, "%d", &ntype);
    printf("number of types  : %d\n", ntype );
    for (int it=0;it<ntype;it++) {
        fscanf(fp,"%s",buffer);
        strcpy(atypei[it],buffer);
        printf("%d %s\n",it,atypei[it]);
    }
    for (int it=0;it<ntype;it++) {
        fscanf(fp,"%d",&data);
        natmi[it]=data;
        printf("%d %d\n",it,natmi[it]);
    }
    fscanf(fp, "%s", buffer);
    for (int i=0;i<MAX_LEN+1;i++){
        cpos[i]=buffer[i];
    }
    printf("%s\n",buffer);
    for (int ia=0;ia<nion;ia++) {
        fscanf(fp,"%s %lf %lf %lf",buffer,&f1,&f2,&f3);
        strcpy(atype[ia],buffer);
        rx[ia]=f1;ry[ia]=f2;rz[ia]=f3;
//       printf("reading %d %s %lf %lf %lf\n",ia,atype[ia],rx[ia],ry[ia],rz[ia]);
    }
    if (strcmp(cpos,"Direct") == 0 ) {
        printf("Direct to Cartesian\n");
        dirkar(nion, rx, ry, rz , Cell->A);
    }
    if (fclose(fp))     { 
       printf("error closing file."); 
       exit(-1); 
    }

    
//    for(int ia=0; ia<nion; ia++) {
//        printf("(after dirkar) reading %d %s %lf %lf %lf\n",ia,atype[ia],rx[ia],ry[ia],rz[ia]);
//    }


    printf("leaving read_posfff\n");
    return 0;
}
