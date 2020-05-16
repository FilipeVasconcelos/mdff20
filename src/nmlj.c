#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "color_mdff.h"
#include "math_mdff.h"
#include "config.h"
#include "field.h"
#include "cell.h"
#include "nmlj.h"
#include "md.h"
#include "io.h"


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
    *rxij = sxij * psimu_cell->A[0][0] + syij * psimu_cell->A[1][0] + szij * psimu_cell->A[2][0];
    *ryij = sxij * psimu_cell->A[0][1] + syij * psimu_cell->A[1][1] + szij * psimu_cell->A[2][1];
    *rzij = sxij * psimu_cell->A[0][2] + syij * psimu_cell->A[1][2] + szij * psimu_cell->A[2][2];

}

int read_nmlj(char* controlfn){

    char buffer[MAX_LEN+1];
    double data[ntype][ntype];
    FILE * fp;
    fp = fopen (controlfn, "r");
    if (NULL == fp )  {
       perror("opening database file");
       return (-1);
    }
    while (EOF != fscanf(fp, "%s", buffer)) {
        // SIGMALJ
        if (strcmp(buffer,"sigmalj") == 0 ) {
            for(int it=0;it<ntype;it++){
                for(int jt=0;jt<ntype;jt++){
                    if (jt >= it) {
                        fscanf(fp,"%lf",&data[it][jt]);
                        sigmalj[it][jt]=data[it][jt];
                    }
                }
            }
        }
        // EPSLJ 
        if (strcmp(buffer,"epslj") == 0 ) {
            for(int it=0;it<ntype;it++){
                for(int jt=0;jt<ntype;jt++){
                    if (jt >= it) {
                        fscanf(fp,"%lf",&data[it][jt]);
                        epslj[it][jt]=data[it][jt];
                    }
                }
            }
        } 
        // PLJ 
        if (strcmp(buffer,"plj") == 0 ) {
            for(int it=0;it<ntype;it++){
                for(int jt=0;jt<ntype;jt++){
                    if (jt >= it) {
                        fscanf(fp,"%lf",&data[it][jt]);
                        plj[it][jt]=data[it][jt];
                    }
                }
            }
        } 
        // QLJ 
        if (strcmp(buffer,"qlj") == 0 ) {
            for(int it=0;it<ntype;it++){
                for(int jt=0;jt<ntype;jt++){
                    if (jt >= it) {
                        fscanf(fp,"%lf",&data[it][jt]);
                        qlj[it][jt]=data[it][jt];
                    }
                }
            }
        } 
        //  lsymmetric 
        if (strcmp(buffer,"lsymmetric") == 0 ) {
            lsymmetric=true;
        } 
        //  trunc 
        if (strcmp(buffer,"trunc") == 0 ) {
            fscanf(fp,"%s",buffer);
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


void init_nmlj(char* controlfn){

    double rcut=cutshortrange;
    rcutsq=rcut*rcut;
    double ppqq[ntype][ntype];
    double one16 = (1.0 / 6.0) ;
    double two16 = pow(2.0,one16);  
    strcpy(trunclabel[0],"no");
    strcpy(trunclabel[1],"linear");
    strcpy(trunclabel[2],"quadratic");

    default_nmlj();

    //read parameters
    read_nmlj(controlfn) ;

    for ( int it=0;it<ntype; it++) {
        for ( int jt=0; jt<ntype; jt++) {
            // symmetric pot
            if (jt>=it){
                plj[jt][it]=plj[it][jt];
                qlj[jt][it]=qlj[it][jt];
                sigmalj[jt][it]=sigmalj[it][jt];
                epslj[jt][it]=epslj[it][jt];
            }
            // p*q
            ppqq[it][jt] = plj[it][jt] * qlj[it][jt];
            // eps / q - p
            epsp[it][jt] = epslj[it][jt] /( qlj [it][jt]- plj[it][jt] );
            // sigma*^2
            sigsq[it][jt]= two16 * two16 * sigmalj[it][jt] * sigmalj[it][jt] ;
            // for the virial                                                          
            fc[it][jt] =  (ppqq [it][jt] * epsp [it][jt]) /  sigsq [it][jt];
            ptwo[it][jt]=plj[it][jt]*0.5;
            qtwo[it][jt]=qlj[it][jt]*0.5;

        }
    }

    info_nmlj();

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
    for(int ia=0;ia<nion;ia++){
        fx[ia]=0.0;fy[ia]=0.0;fz[ia]=0.0;
    }
    for (int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            tau_nonb[i][j]=0.0;
        }
    }
    /*************************************** 
            cartesian to direct                    
     ***************************************/
//    printf("before kardir %f %f\n",ry[100],psimu_cell->B[0][0]);
    kardir ( nion , rx , ry , rz , psimu_cell->B ) ;
//    printf("after kardir %f %f\n",ry[100],psimu_cell->B[0][0]);

//    printf("inside engforce_nmlj_pbc %d\n",nion);
    
    //printf("%d %d\n",patom_dec->iaStart,patom_dec->iaEnd);
    for(int ia=patom_dec->iaStart;ia<patom_dec->iaEnd;ia++) {
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
		    fx[ia] += fxij;
                    fy[ia] += fyij;
                    fz[ia] += fzij;
                    fx[ja] -= fxij;
                    fy[ja] -= fyij;
                    fz[ja] -= fzij;
                    // stress tensor
                    tau_nonb[0][0] += (rxij*fxij + rxij*fxij)*0.5;
                    tau_nonb[0][1] += (rxij*fyij + ryij*fxij)*0.5;
                    tau_nonb[0][2] += (rxij*fzij + rzij*fxij)*0.5;
                    tau_nonb[1][0] += (ryij*fxij + rxij*fyij)*0.5;
                    tau_nonb[1][1] += (ryij*fyij + ryij*fyij)*0.5;
                    tau_nonb[1][2] += (ryij*fzij + rzij*fyij)*0.5;
                    tau_nonb[2][0] += (rzij*fxij + rxij*fzij)*0.5;
                    tau_nonb[2][1] += (rzij*fyij + ryij*fzij)*0.5;
                    tau_nonb[2][2] += (rzij*fzij + rzij*fzij)*0.5;
                }
            }
        }
    }
    double coeff=1./(psimu_cell->Omega*press_unit);
    for (int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            tau_nonb[i][j] *= coeff;
        }
    }

    /*************************************** 
            direct to cartesian                   
     ***************************************/
    dirkar ( nion , rx , ry , rz , psimu_cell->A ) ;

    //if (itime%nprint==0) printf("u = %f\n",*u);

}

void info_nmlj(){

    if (ionode) {
        LSEPARATOR;
        printf("n-m lennard-jones\n");
        LSEPARATOR;
        putchar('\n');
        printf("       eps    /    / sigma* \\ q       / sigma*  \\ p  \\ \n");
        printf(" V = ------- |  p | ------- |    - q | -------- |    |\n"); 
        printf("      q - p   \\    \\   r    /         \\    r    /    /\n"); 
    
        putchar('\n');
        putchar('\n');
        printf("cutoff      = %12.4f \n",cutshortrange);
        printf("truncation  = %s (%d)\n",trunclabel[trunctype],trunctype);        
        putchar('\n');
        printf("pair interactions\n");
        for(int it=0;it<ntype;it++){
            for(int jt=0;jt<ntype;jt++){
                if( (jt>=it) && (jt==it)){
                    LSEPARATOR;
                    printf(" %s <--> %s  \n",atypei[it],atypei[jt]);
                    LSEPARATOR;
                    printf("sigma = %16.8e\n",sigmalj[it][jt]);
                    printf("eps   = %16.8e\n",epslj[it][jt]);
                    printf("q     = %16.8e\n",qlj[it][jt]);
                    printf("p     = %16.8e\n",plj[it][jt]);
                }
                else if( (jt>it) ){
                    LSEPARATOR;
                    printf(" %s <--> %s  \n",atypei[it],atypei[jt]);
                    LSEPARATOR;
                    printf("sigma = %16.8e (%-16.8e)\n",sigmalj[it][jt],sigmalj[jt][it]);
                    printf("eps   = %16.8e (%-16.8e)\n",epslj[it][jt],epslj[jt][it]);
                    printf("q     = %16.8e (%-16.8e)\n",qlj[it][jt],qlj[jt][it]);
                    printf("p     = %16.8e (%-16.8e)\n",plj[it][jt],plj[jt][it]);
                }
            }
        }
        putchar('\n');
    }
}
