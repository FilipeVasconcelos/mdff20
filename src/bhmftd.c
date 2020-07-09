#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "tt_damp.h"
#include "global.h"
#include "field.h"
#include "config.h"
#include "pbc.h"
#include "io.h"
#include "nonbonded.h"
#include "cell.h"
#include "verlet.h"
#include "bhmftd.h"
#include "timing.h"

#ifdef DEBUG
    #define DEBUG_BHMFTD
#endif
//#define DEBUG_
#ifdef DEBUG_
    #define DEBUG_BHMFTD
#endif

/******************************************************************************/
int read_bhmftd(char* controlfn) {

    char buffer[MAX_LEN+1];
    FILE *fp;
    fp = fopen (controlfn, "r");
    if (NULL == fp ) {
       pError("opening control file (reading bhmftd)");
       return (-1);
    }
    while (EOF != fscanf(fp, "%s", buffer)) {

        // Abhmftd
        if (strcmp(buffer,"Abhmftd") == 0 ) {
            for(int it=0;it<ntype;it++){
                for(int jt=0;jt<ntype;jt++){
                    if (jt >= it) {
                        fscanf(fp,"%lf",&Abhmftd[it][jt]);
                    }
                }
            }
        }
        // Bbhmftd
        if (strcmp(buffer,"Bbhmftd") == 0 ) {
            for(int it=0;it<ntype;it++){
                for(int jt=0;jt<ntype;jt++){
                    if (jt >= it) {
                        fscanf(fp,"%lf",&Bbhmftd[it][jt]);
                    }
                }
            }
        }
        // Cbhmftd
        if (strcmp(buffer,"Cbhmftd") == 0 ) {
            for(int it=0;it<ntype;it++){
                for(int jt=0;jt<ntype;jt++){
                    if (jt >= it) {
                        fscanf(fp,"%lf",&Cbhmftd[it][jt]);
                    }
                }
            }
        }
        // Dbhmftd
        if (strcmp(buffer,"Dbhmftd") == 0 ) {
            for(int it=0;it<ntype;it++){
                for(int jt=0;jt<ntype;jt++){
                    if (jt >= it) {
                        fscanf(fp,"%lf",&Dbhmftd[it][jt]);
                    }
                }
            }
        }
        // BDbhmftd
        if (strcmp(buffer,"BDbhmftd") == 0 ) {
            for(int it=0;it<ntype;it++){
                for(int jt=0;jt<ntype;jt++){
                    if (jt >= it) {
                        fscanf(fp,"%lf",&BDbhmftd[it][jt]);
                    }
                }
            }
        }
    }
    return 0;
}

/******************************************************************************/
void init_bhmftd(char* controlfn) {

    //default_bhmftd();
    /* read parameters */
    read_bhmftd(controlfn);

    for ( int it=0;it<ntype; it++) {
        for ( int jt=0; jt<ntype; jt++) {
            // symmetric pot
            if (jt>=it){
                Abhmftd[jt][it] =Abhmftd[it][jt];
                Bbhmftd[jt][it] =Bbhmftd[it][jt];
                Cbhmftd[jt][it] =Cbhmftd[it][jt];
                Dbhmftd[jt][it] =Dbhmftd[it][jt];
                BDbhmftd[jt][it]=BDbhmftd[it][jt];
            }
        }
    }

    if (lbhmftd) get_TT_damp();

    info_bhmftd();
}
/******************************************************************************/
void engforce_bhmftd_pbc(double *u, double *pvir, double tau[3][3])
{
    double rxi,ryi,rzi;
    double rxij,ryij,rzij;
    double fxij,fyij,fzij;
    double rijsq;
    double wij;

    int p1,p2;
    double uu =0.0;
    double ttau[3][3];
    int ia,j1,ja,jb,je;
    double ir2,d,erh,ir6,ir8,ir6d,ir8d,f6,f8,fdiff6,fdiff8,ir7,ir9;

    for (int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            ttau[i][j]=0.0;
        }
    }
    /***************************************
            cartesian to direct
     ***************************************/
    kardir ( nion , rx , ry , rz , simuCell.B ) ;

#ifdef DEBUG_BHMFTD
    int counttest=0;
#endif

    #pragma omp parallel default(none)\
                         shared(srcutsq,rx,ry,rz,vx,vy,vz,fx,fy,fz,ttau,uu,typia,nion,atomDec,lverletL,verlet_nb,\
                                 Abhmftd,Bbhmftd,Cbhmftd,Dbhmftd,BDbhmftd,lbhmftd)\
                         private(ia,j1,jb,je,ja,rxi,ryi,rzi,rxij,ryij,rzij,rijsq,p1,p2,wij,\
                                 fxij,fyij,fzij,ir2,d,erh,ir6,ir8,ir6d,ir8d,f6,f8,fdiff6,fdiff8,ir7,ir9)
    {
        #pragma omp for reduction (+:uu,ttau,fx[:nion],fy[:nion],fz[:nion]) schedule(dynamic,16)
        for(ia=atomDec.iaStart;ia<atomDec.iaEnd;ia++) {
            rxi = rx[ia];
            ryi = ry[ia];
            rzi = rz[ia];
            if (lverletL) {
                jb=verlet_nb->point[ia];
                je=verlet_nb->point[ia+1];
            }
            else{
                jb = 0;
                je = nion;
            }
            for (j1=jb;j1<je;j1++) {
                if (lverletL){
                    ja=verlet_nb->list[j1];
                }
                else{
                    ja=j1;
                }
                if ( (( ja > ia ) && !lverletL ) ||
                     (( ja !=ia ) && lverletL) )  {
                    p1 = typia[ia];
                    p2 = typia[ja];
                    if ( Abhmftd[p1][p2] == 0.0 ) continue;
                    rxij = rxi - rx[ja];
                    ryij = ryi - ry[ja];
                    rzij = rzi - rz[ja];
                    pbc(&rxij,&ryij,&rzij);
                    rijsq = rxij * rxij + ryij * ryij + rzij * rzij;

                    if ( rijsq < srcutsq ) {
                        ir2 = 1.0 / rijsq;
                        d = sqrt(rijsq);
                        erh = Abhmftd[p1][p2] * exp( - Bbhmftd[p1][p2] * d ) ;
                        ir6 = ir2 * ir2 * ir2;
                        ir8 = ir6 * ir2;
                        ir6 = ir6 * Cbhmftd[p1][p2];
                        ir8 = ir8 * Dbhmftd[p1][p2];
                        /* damping TT */
                        if (lbhmftd) {
                            TT_damping_functions(BDbhmftd[p1][p2],1.0 ,d ,&f6 ,&fdiff6, 6);
                            TT_damping_functions(BDbhmftd[p1][p2],1.0 ,d ,&f8 ,&fdiff8, 8);
                            ir6d=ir6*f6;
                            ir8d=ir8*f8;
                        }
                        else
                        {
                            fdiff6 = 0.0;
                            fdiff8 = 0.0;
                            ir6d=ir6;
                            ir8d=ir8;
                        }
                        ir7 = 6.0 * ir6d / d;
                        ir9 = 8.0 * ir8d / d;
#ifdef DEBUG_BHMFTD
            counttest+=1;
            io_node printf("in bhmft main loop %d %d %e %e %e %e\n",ia,ja,uu,erh - ir6d - ir8d,sqrt(rijsq),Abhmftd[p1][p2]);
#endif
                        uu += erh - ir6d - ir8d;
                        wij  =  Bbhmftd[p1][p2] * erh - ir7 - ir9 + ir6 * fdiff6 + ir8 * fdiff8;
                        wij  = wij / d;
                        fxij = wij * rxij;
                        fyij = wij * ryij;
                        fzij = wij * rzij;
                        // forces
                        fx[ia] +=  fxij;
                        fy[ia] +=  fyij;
                        fz[ia] +=  fzij;
                        fx[ja] += -fxij;
                        fy[ja] += -fyij;
                        fz[ja] += -fzij;
                        // stress tensor (symmetric!)
                        ttau[0][0] += (rxij*fxij + rxij*fxij)*0.5;
                        ttau[0][1] += (rxij*fyij + ryij*fxij)*0.5;
                        ttau[0][2] += (rxij*fzij + rzij*fxij)*0.5;
                        ttau[1][1] += (ryij*fyij + ryij*fyij)*0.5;
                        ttau[1][2] += (ryij*fzij + rzij*fyij)*0.5;
                        ttau[2][2] += (rzij*fzij + rzij*fzij)*0.5;
                    }
                }
            }
        }
    }
    ttau[1][0]=ttau[0][1];
    ttau[2][0]=ttau[0][2];
    ttau[2][1]=ttau[1][2];
    for (int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            tau[i][j]=ttau[i][j]*simuCell.inveOmegaPU;
        }
    }
    /* Potential energy */
    *u=uu;
    /* Pvirial */
    *pvir=0.0;
    for (int i=0;i<3;i++){
        *pvir+=tau[i][i]*onethird;
    }
#ifdef MPI
    statime(8);
    MPI_Allreduce_sumDouble(u,1);
    MPI_Allreduce_sumDouble(pvir,1);
    MPI_Allreduce_sumDouble(fx,nion);
    MPI_Allreduce_sumDouble(fy,nion);
    MPI_Allreduce_sumDouble(fz,nion);
    MPI_Allreduce_sumDouble(tau[0],3);
    MPI_Allreduce_sumDouble(tau[1],3);
    MPI_Allreduce_sumDouble(tau[2],3);
    statime(9);
    mestime(&COMMCPUtime,9,8);
#endif
    /***************************************
            direct to cartesian
     ***************************************/
    dirkar ( nion , rx , ry , rz , simuCell.A ) ;
#ifdef DEBUG_BHMFTD
    io_node printf("exiting bhmft engforce count pairs %d %d\n",counttest,(nion*(nion-1))/2);
#endif

}
/******************************************************************************/
void info_bhmftd(){
    if (ionode) {
        SEPARATOR;
        printf("BHMFT field info\n");
        LSEPARATOR;
        if (lbhmftd) {
            printf("Born-Huggins-Meyer-Fumi-Tossi + Tang-Toennies Damping \n");
            LSEPARATOR;
            putchar('\n');
            printf("                                  C              D  \n");
            printf(" V = A  exp ( - B  r ) - f_6(r) ----- - f_8(r) -----\n");
            printf("                                 r^6            r^8 \n");
            putchar('\n');
            putchar('\n');
            printf("                             n                      \n");
            printf("                            ----   (BD * r)^k       \n");
            printf(" f_n(r) = 1 - exp(-BD r ) * \\     ------------     \n");
            printf("                            /___       k!           \n");
            printf("                             k=0                    \n");
        }
        else
        {
            printf("Born-Huggins-Meyer-Fumi-Tossi \n");
            LSEPARATOR;
            putchar('\n');
            printf("                            C        D  \n");
            printf(" V = A  exp ( - B  r ) -  ----- -  -----\n");
            printf("                           r^6      r^8 \n");
        }

        putchar('\n');
        putchar('\n');
        printf("cutoff      = %-6.2f \n",cutshortrange);
//        printf("truncation  = %s (%d)\n",trunclabel[trunctype],trunctype);
//        printf("long range correction (energy)   : "ee"\n",utail);
//        printf("long range correction (pressure) : "ee"\n",ptail);
        LSEPARATOR;
        printf("pair interactions\n");
        for(int it=0;it<ntype;it++){
            for(int jt=0;jt<ntype;jt++){
                if( (jt>=it) && (jt==it)){
                    LSEPARATOR;
                    printf(" %s <--> %s  \n",atypit[it],atypit[jt]);
                    LSEPARATOR;
                    printf("A      = %16.8e\n",Abhmftd[it][jt]);
                    printf("B      = %16.8e\n",Bbhmftd[it][jt]);
                    printf("C      = %16.8e\n",Cbhmftd[it][jt]);
                    printf("D      = %16.8e\n",Dbhmftd[it][jt]);
                    if (lbhmftd) printf("BD     = %16.8e\n",BDbhmftd[it][jt]);
                }
                else if( (jt>it) ){
                    LSEPARATOR;
                    printf(" %s <--> %s  \n",atypit[it],atypit[jt]);
                    LSEPARATOR;
                    printf("A      = %16.8e (%-16.8e)\n",Abhmftd[it][jt],Abhmftd[jt][it]);
                    printf("B      = %16.8e (%-16.8e)\n",Bbhmftd[it][jt],Bbhmftd[jt][it]);
                    printf("C      = %16.8e (%-16.8e)\n",Cbhmftd[it][jt],Cbhmftd[jt][it]);
                    printf("D      = %16.8e (%-16.8e)\n",Dbhmftd[it][jt],Dbhmftd[jt][it]);
                    if (lbhmftd) printf("BD     = %16.8e (%-16.8e)\n",BDbhmftd[it][jt],BDbhmftd[jt][it]);
                }
            }
        }
        putchar('\n');
    } /* ionode */
}

