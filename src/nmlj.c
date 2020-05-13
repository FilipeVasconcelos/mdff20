#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "config.h"
#include "field.h"
#include "cell.h"
#include "nmlj.h"
#include "md.h"

//#define nint(x) ( (x) >=0  ? (int)((x)+0.5) : (int) ((x)-0.5) )
//#define nint(x) ( (x) >= 0.5  ? (int) ((x)+0.5) : (int) ((x)-0.5) )
#define nint(x) nearbyint(x)

double addtruncU(int p1 , int p2, double srp, double srq, int trunc){
    if (trunc == 0) {
        return epsp[p1][p2] * ( plj[p1][p2] * srq - qlj[p1][p2] * srp ) ;
    }
    return 0.0;
}

void pbc(double *rxij,double *ryij,double *rzij){
    double sxij = *rxij - nint ( *rxij );
    double syij = *ryij - nint ( *ryij );
    double szij = *rzij - nint ( *rzij );
    //printf("test nint %f %f %f %f %f %f %d %d %d",sxij,syij,szij,*rxij,*ryij,*rzij,nint ( *rxij ),nint ( *ryij ),nint ( *rzij ));
    *rxij = sxij * psimu_cell->A[0][0] + syij * psimu_cell->A[0][1] + szij * psimu_cell->A[0][2];
    *ryij = sxij * psimu_cell->A[1][0] + syij * psimu_cell->A[1][1] + szij * psimu_cell->A[1][2];
    *rzij = sxij * psimu_cell->A[2][0] + syij * psimu_cell->A[2][1] + szij * psimu_cell->A[2][2];
    //printf(" after : %f %f %f %f\n",*rxij,*ryij,*rzij,psimu_cell->A[0][0]);

}

void read_nmlj(){


}

void init_nmlj(){
    double rcut=cutshortrange;
    rcutsq=rcut*rcut;
    double rcut3=rcutsq*rcut;
    double rskin = rcut + skindiff;
    double rskinsq = rskin*rskin;
    
    double ppqq[ntype][ntype];
    double one13 = (1.0 / 3.0) ;
    double one16 = (1.0 / 6.0) ;
    double two16 = pow(2.0,one16);  

/*
 * Temp 
*/
    for ( int it=0;it<ntype; it++){
        for ( int jt=0; jt<ntype; jt++){
            printf("pair potential %d <--> %d\n",it,jt);
            plj[it][jt]=6.0;
            qlj[it][jt]=12.0;
            ptwo[it][jt]=plj[it][jt]*0.5;
            qtwo[it][jt]=qlj[it][jt]*0.5;
            sigmalj[it][jt]=3.405;
            epslj[it][jt]=0.010323576;
            //sigmalj[it][jt]=1.0;
            //epslj[it][jt]=1.0;
        }
    }
//
//
    for ( int it=0;it<ntype; it++)    {
        for ( int jt=0; jt<ntype; jt++)        {
            printf("pair potential %d <--> %d",it,jt);
            ppqq[it][jt] = plj[it][jt] * qlj[it][jt];
            // eps / q - p
            epsp[it][jt] = epslj[it][jt] /( qlj [it][jt]- plj[it][jt] );
            // sigma*^2
            sigsq[it][jt]= two16 * two16 * sigmalj[it][jt] * sigmalj[it][jt] ;
            // for the virial                                                          
            fc[it][jt] =  ppqq [it][jt] * epsp [it][jt] /  sigsq [it][jt];
            printf("fc : %f sigsq : %f epsp : %f pq : %f\n",fc[it][jt],sigsq[it][jt],epsp[it][jt],ppqq[it][jt]);

        }
    }
}

void engforce_nmlj_pbc(double *u)
{
    double rxi,ryi,rzi;
    double rxij,ryij,rzij;
    double fxij,fyij,fzij;
    double rijsq;
    double wij,vir;
    double sr2,srp,srq;

    int p1,p2;
    *u=0;

    /*************************************** 
            cartesian to direct                    
     ***************************************/
//    printf("before kardir %f %f\n",ry[100],psimu_cell->B[0][0]);
    kardir ( nion , rx , ry , rz , psimu_cell->B ) ;
//    printf("after kardir %f %f\n",ry[100],psimu_cell->B[0][0]);

//    printf("inside engforce_nmlj_pbc %d\n",nion);

    for(int ia=0;ia<nion;ia++) {
        rxi = rx[ia];
        ryi = ry[ia];
        rzi = rz[ia];

        int jb = 0;
        int je = nion;
        for (int ja=jb;ja<je;ja++) {
            if ( ja > ia ) {
                rxij = rxi - rx[ja];
                ryij = ryi - ry[ja];
                rzij = rzi - rz[ja];
                pbc(&rxij,&ryij,&rzij);
    //            printf("%f %f %f\n",rxij,ryij,rzij);      
                rijsq = rxij * rxij + ryij * ryij + rzij * rzij;
                p1 = itype[ia];
                p2 = itype[ja];
    
                    //printf("%d %d %f %f %f %f %d %d %f %f %f\n",ia,ja,rijsq,rcutsq,cutshortrange,*u,p1,p2,srp,srq,sr2);
                if ( rijsq < rcutsq ) {
                    sr2 = sigsq[p1][p2] / rijsq;
                    srp = pow(sr2,ptwo[p1][p2]);
                    srq = pow(sr2,qtwo[p1][p2]);
                    *u+=addtruncU(p1,p2,srp,srq,0);
    //                printf("%d %d %f %f %f %f %d %d %f %f %f %f %f\n",ia,ja,rijsq,rcutsq,cutshortrange,*u,p1,p2,srp,srq,sr2,ptwo[p1][p2],qtwo[p1][p2]);
    //                if (ia == 2) exit(0);
		
                    wij = fc[p1][p2] * (srq-srp) * sr2;
                    fxij = wij * rxij;
                    fyij = wij * ryij;
                    fzij = wij * rzij;
                    // virial;
                    vir = vir + wij * rijsq;
                    // forces
		    fx [ ia ] = fx [ ia ] + fxij;
                    fy [ ia ] = fy [ ia ] + fyij;
                    fz [ ia ] = fz [ ia ] + fzij;
                    fx [ ja ] = fx [ ja ] - fxij;
                    fy [ ja ] = fy [ ja ] - fyij;
                    fz [ ja ] = fz [ ja ] - fzij;
                }
            }
        }
    }

    /*************************************** 
            direct to cartesian                   
     ***************************************/
    dirkar ( nion , rx , ry , rz , psimu_cell->A ) ;

    //if (itime%nprint==0) printf("u = %f\n",*u);

}

