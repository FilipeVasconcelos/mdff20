#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "color_mdff.h"
#include "config.h"
#include "field.h"
#include "timing.h"
#include "cell.h"
#include "nmlj.h"
#include "md.h"
#include "io.h"
#include "pbc.h"
#include "global.h"
#include "verlet.h"
#include <omp.h>

#ifdef DEBUG
    #define DEBUG_NMLJ
#endif
//#define DEBUG_
#ifdef DEBUG_
    #define DEBUG_NMLJ
#endif
/******************************************************************************/
double addtruncU(int p1 , int p2, double srp, double srq){
    if (trunctype == 0) {
        return epsqp[p1][p2] * ( plj[p1][p2] * srq - qlj[p1][p2] * srp ) ;
    }
    if (trunctype == 1) {
        return epsqp[p1][p2] * ( plj[p1][p2] * srq - qlj[p1][p2] * srp ) - uc[p1][p2];
    }
    return 0.0;
}

/******************************************************************************/
int read_nmlj(char* controlfn){

    int c;
    char buffer[MAX_LEN+1];
    FILE * fp;
    fp = fopen (controlfn, "r");
    if (NULL == fp )  {
        pError("opening control file (reading nmlj)");
        return (-1);
    }
    while (EOF != fscanf(fp, "%s", buffer)) {
        // SIGMALJ
        if (strcmp(buffer,"sigmalj") == 0 ) {
            for(int it=0;it<ntype;it++){
                for(int jt=0;jt<ntype;jt++){
                    if (jt >= it) {
                        c=fscanf(fp,"%lf",&sigmalj[it][jt]);
                    }
                }
            }
        }
        // EPSLJ
        if (strcmp(buffer,"epslj") == 0 ) {
            for(int it=0;it<ntype;it++){
                for(int jt=0;jt<ntype;jt++){
                    if (jt >= it) {
                        c=fscanf(fp,"%lf",&epslj[it][jt]);
                    }
                }
            }
        }
        // PLJ
        if (strcmp(buffer,"plj") == 0 ) {
            for(int it=0;it<ntype;it++){
                for(int jt=0;jt<ntype;jt++){
                    if (jt >= it) {
                        c=fscanf(fp,"%lf",&plj[it][jt]);
                    }
                }
            }
        }
        // QLJ
        if (strcmp(buffer,"qlj") == 0 ) {
            for(int it=0;it<ntype;it++){
                for(int jt=0;jt<ntype;jt++){
                    if (jt >= it) {
                        c=fscanf(fp,"%lf",&qlj[it][jt]);
                    }
                }
            }
        }
        //  trunc
        if (strcmp(buffer,"trunc") == 0 ) {
            c=fscanf(fp,"%s",buffer);
            trunctype=-1;
            for (int k=0;k<3;k++){
                if (strcmp(buffer,trunclabel[k]) == 0){
                    trunctype=k;
                }
            }
            if (trunctype<0){
                pError("trunc not allowed\n");
                printf("%s\n",buffer);
                exit(-1);
            }
        }
    }
    return 0;
}

/******************************************************************************/
void default_nmlj(){
    for ( int it=0;it<ntype; it++){
        for ( int jt=0; jt<ntype; jt++){
            plj[it][jt]=6.0;
            qlj[it][jt]=12.0;
            sigmalj[it][jt]=3.405;
            epslj[it][jt]=.010323576;
        }
    }
}


/******************************************************************************/
void init_nmlj(char* controlfn){

    double srcut3,pp,qq,ppqq,ut,pt,epsppqq,srq,srp,sr,sr2,qq3,pp3;
    double one16 = (1.0 / 6.0) ;
    double two16 = pow(2.0,one16);
    strcpy(trunclabel[0],"no");
    strcpy(trunclabel[1],"linear");
    strcpy(trunclabel[2],"quadratic");

    default_nmlj();

    //read parameters
    read_nmlj(controlfn) ;

/*
 ==================================================================================================
 TAIL ( checked september 2011) :
 be careful two minus are vanishing

     2 pi rc^3 espilon NA NB    /     p         / sigma* \ q           q         / sigma*  \ p  \
    -------------------------- |  ----------   | ------- |    --   ----------   | -------- |    |
         ( q - p )   V          \ ( q - 3 )     \   rc   /          ( p - 3 )    \  rc     /    /

  virial the same with qp in the first term

 ==================================================================================================
*/

/*
 ==================================================================================================
  simple truncation

  trunclabel = 'linear'
  trunctype  = 1


           eps    /    / sigma* \ q         / sigma* \ p  \
  V  =   ------- |  p | ------- |    -  q  | --------|    |   -  c1   with  sigma* = 2^(1/6)*sigma
          q - p   \    \   r    /           \    r   /    /


          eps      /     /  sigma* \ q        /  sigma* \ p  \
  c1 = ---------  |   p |  -------- |    - q |  -------- |   |      with rc = cutshortrange
         q - p     \     \   rc    /          \   rc    /    /

 ==================================================================================================
*/

/*
 ==================================================================================================
  truncation presented in J. Chem. Phys. 135 (2011) , Sengupta, Vasconcelos, Affouard, Sastry

  trunclabel = 'quadratic'
  trunctype  = 2


         eps    /    / sigma* \ q         / sigma* \ p  \
V  =   ------- |  p | ------- |    -  q  | --------|    |   +  c1 r^2  -  c2 with sigma* 2^(1/6)*sigma
        q - p   \    \   r   /            \    r   /    /


 verified on Sage (January 2013)

              eps p  q           /  / sigma* \ q       /  sigma* \ p \
    c1 =  --------------------- |  | --------|    -   | --------- |  |
          2  ( q - p ) * rc^2    \  \   rc   /         \   rc    /   /

   and

           epsilon     /           / sigma* \ q               / sigma* \ p  \
    c2 =  ----------- |  (2p+qp)  | --------|     - (2q+pq)  | --------|   |
          2 ( q - p )  \           \   rc   /                 \   rc   /   /

 ==================================================================================================
 */
    for ( int it=0;it<ntype; it++) {
        for ( int jt=0; jt<ntype; jt++) {
            // symmetric pot
            if (jt>=it){
                plj[jt][it]=plj[it][jt];
                qlj[jt][it]=qlj[it][jt];
                sigmalj[jt][it]=sigmalj[it][jt];
                epslj[jt][it]=epslj[it][jt];
            }
            srcut3=srcutsq*cutshortrange;
            pp=plj[it][jt];qq=qlj[it][jt];
            pp3=pp-3.0; qq3=qq-3.0;

            /* -------------------------------------------------------------- */
            // used in main energy force function
            // eps / q - p
            epsqp[it][jt] = epslj[it][jt] /(qq - pp);
            // sigma*^2
            sigsq[it][jt]= two16 * two16 * sigmalj[it][jt] * sigmalj[it][jt] ;
            // for the virial
            fc[it][jt] =  (pp * qq * epsqp [it][jt]) / sigsq [it][jt];
            ptwo[it][jt]=pp*0.5;
            qtwo[it][jt]=qq*0.5;
            /* -------------------------------------------------------------- */

            // local
            // p*q
            epsppqq=epsqp[it][jt];
            ppqq = pp * qq;
            sr2=sigsq[it][jt]/srcutsq;
            sr=sqrt(sr2);
            srp=pow(sr,pp);
            srq=pow(sr,qq);
            //trunctype = 1
            uc[it][jt]=epsppqq*(pp*srq-qq*srp);
            //trunctype = 2
            c1[it][jt]=epsppqq*ppqq/(2.0*srcutsq) * (srq-srp);
            c2[it][jt]=0.5*epsppqq*((2.0*pp+ppqq)*srq-
                                 (2.0*qq+ppqq)*srp);
            // tail energy
            ut=epsppqq*(pp*srq/qq3 - qq*srp/pp3);
            pt=ppqq*epsppqq*(srq/qq3 -srp/pp3);
            ut*=srcut3*TPI;
            pt*=srcut3*TPI;
            if ( (nionit[it] != 0) && (nionit[jt] != 0) ) {
                utail+=ut*nionit[it]*nionit[jt]*simuCell.inveOmega;
                ptail+=pt*nionit[it]*nionit[jt]*simuCell.inveOmega;
            }
        }
    }
    ptail=ptail/(press_unit*3.0); /* /inveOmega in mdff(fortran) */

    info_nmlj();
}

/******************************************************************************/
void engforce_nmlj_pbc(double *u, double *pvir, double tau[3][3])
{
    double rxi,ryi,rzi;
    double rxij,ryij,rzij;
    double fxij,fyij,fzij;
    double rijsq;
    double wij;
    double sr2,srp,srq;

    int p1,p2;
    double uu =0.0;
//    double pvv =0.0;
    double ttau[3][3];
    int ia,j1,ja,jb,je;

    for (int i=0;i<3;i++){
        for(int j=i;j<3;j++){
            ttau[i][j]=0.0;
        }
    }
    /***************************************
            cartesian to direct
     ***************************************/
    kardir ( nion , rx , ry , rz , simuCell.B ) ;

#ifdef DEBUG_NMLJ
    int counttest=0;
#endif

    #pragma omp parallel default(none) \
                         shared(srcutsq,rx,ry,rz,vx,vy,vz,fx,fy,fz,sigsq,ptwo,qtwo,typia,epsqp,plj,qlj,fc,\
                                uc,nion,ttau,uu,atomDec,lverletL,verlet_nb) \
                         private(ia,j1,jb,je,ja,rxi,ryi,rzi,rxij,ryij,rzij,rijsq,p1,p2,sr2,srp,srq,wij,fxij,fyij,fzij)
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
            else {
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
#ifdef DEBUG_NMLJ
                    counttest+=1;
                    io_node printf("in nmlj main loop %d %d %d %d\n",ia,ja,jb,je);
#endif
                    rxij = rxi - rx[ja];
                    ryij = ryi - ry[ja];
                    rzij = rzi - rz[ja];
                    pbc(&rxij,&ryij,&rzij);
                    rijsq = rxij * rxij + ryij * ryij + rzij * rzij;
                    if ( rijsq < srcutsq ) {
                        p1 = typia[ia];
                        p2 = typia[ja];
                        sr2 = sigsq[p1][p2] / rijsq;
                        srp = pow(sr2,ptwo[p1][p2]);
                        srq = pow(sr2,qtwo[p1][p2]);
                        uu += addtruncU(p1,p2,srp,srq);
                        wij = fc[p1][p2] * (srq-srp) * sr2;
                        fxij = wij * rxij;
                        fyij = wij * ryij;
                        fzij = wij * rzij;
                        // Pvirial;
                        //pvv += wij * rijsq;
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
    *u=uu;
    *pvir=0.0;
    for (int i=0;i<3;i++){
        *pvir += tau[i][i] * onethird;
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
#ifdef DEBUG_NMLJ
    io_node printf("exiting nmlj engforce count pairs %d %d\n",counttest,(nion*(nion-1))/2);
#endif

}

/******************************************************************************/
void info_nmlj(){

    if (ionode) {
        SEPARATOR;
        printf("nmlj field info\n");
        LSEPARATOR;
        printf("n-m lennard-jones\n");
        LSEPARATOR;
        putchar('\n');
        printf("       eps    /    / sigma* \\ q       / sigma*  \\ p  \\ \n");
        printf(" V = ------- |  p | ------- |    - q | -------- |    |\n");
        printf("      q - p   \\    \\   r    /         \\    r    /    /\n");

        putchar('\n');
        putchar('\n');
        printf("cutoff      = %-6.2f \n",cutshortrange);
        printf("truncation  = %s (%d)\n",trunclabel[trunctype],trunctype);
        printf("long range correction (energy)   : "ee"\n",utail);
        printf("long range correction (pressure) : "ee"\n",ptail);
        LSEPARATOR;
        printf("pair interactions\n");
        for(int it=0;it<ntype;it++){
            for(int jt=0;jt<ntype;jt++){
                if( (jt>=it) && (jt==it)){
                    LSEPARATOR;
                    printf(" %s <--> %s  \n",atypit[it],atypit[jt]);
                    LSEPARATOR;
                    printf("sigma = %16.8e\n",sigmalj[it][jt]);
                    printf("eps   = %16.8e\n",epslj[it][jt]);
                    printf("q     = %16.8e\n",qlj[it][jt]);
                    printf("p     = %16.8e\n",plj[it][jt]);
                }
                else if( (jt>it) ){
                    LSEPARATOR;
                    printf(" %s <--> %s  \n",atypit[it],atypit[jt]);
                    LSEPARATOR;
                    printf("sigma = %16.8e (%-16.8e)\n",sigmalj[it][jt],sigmalj[jt][it]);
                    printf("eps   = %16.8e (%-16.8e)\n",epslj[it][jt],epslj[jt][it]);
                    printf("q     = %16.8e (%-16.8e)\n",qlj[it][jt],qlj[jt][it]);
                    printf("p     = %16.8e (%-16.8e)\n",plj[it][jt],plj[jt][it]);
                }
            }
        }
        putchar('\n');
    } /* ionode */
}
