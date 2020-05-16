#include <stdio.h>
#include <string.h>
#include <math.h>

#include "cell.h"
#include "constants.h"
#include "tools.h"
#include "io.h"

void lattice(CELL * Cell)
{
    Cell->G[0][0] = Cell->A[0][0]*Cell->A[0][0] + 
                    Cell->A[1][0]*Cell->A[1][0] + 
                    Cell->A[2][0]*Cell->A[2][0] ;
    Cell->G[1][0] = Cell->A[1][0]*Cell->A[1][1] + 
                    Cell->A[2][0]*Cell->A[2][1];
    Cell->G[2][0] = Cell->A[2][0]*Cell->A[2][2];
    Cell->G[1][1] = Cell->A[1][1]*Cell->A[1][1] + 
                    Cell->A[2][2]*Cell->A[2][2];
    Cell->G[2][1] = Cell->A[2][1]*Cell->A[2][2];
    Cell->G[0][2] = Cell->G[2][0];
    Cell->G[1][2] = Cell->G[2][1];
    Cell->G[2][2] = Cell->A[2][2]*Cell->A[2][2];

    expro(Cell->B[0],Cell->A[1],Cell->A[2]);
    expro(Cell->B[1],Cell->A[2],Cell->A[0]);
    expro(Cell->B[2],Cell->A[0],Cell->A[1]);
    
    
    // volume ( direct )
    double omega;
    omega = Cell->B[0][0]*Cell->A[0][0] +
            Cell->B[1][0]*Cell->A[1][0] + 
            Cell->B[2][0]*Cell->A[2][0];
    Cell->Omega=omega;
    Cell->inveOmega=1.0/omega;

    // shortest distance between opposite faces                                                                          
    for(int i=0;i<3;i++){
        Cell->w[i] = omega/sqrt(Cell->B[0][i]*Cell->B[0][i] + 
                                Cell->B[1][i]*Cell->B[1][i] + 
                                Cell->B[2][i]*Cell->B[2][i]);
    }

    for(int i=0;i<3;i++) {
        for(int j=0;j<3;j++) {
            Cell->B[i][j]=Cell->B[i][j] * Cell->inveOmega ;
        }
    }
    // norms (direct and reciprocal)
    for (int i=0;i<3;i++) {
        for(int j=0;j<3;j++) {
            Cell->Anorm[i]+=Cell->A[j][i]*Cell->A[j][i];
            Cell->Bnorm[i]+=Cell->B[j][i]*Cell->B[j][i];
        }
        Cell->Anorm[i]=sqrt(Cell->Anorm[i]);
        Cell->Bnorm[i]=sqrt(Cell->Bnorm[i]);
    }
    // angles (direct)
    double ang[3];
    angles_(&ang,Cell->A,Cell->Anorm);
    for(int i=0;i<3;i++){
        Cell->ang[i]=ang[i];
    }
    // volume ( reciprocal )
    Cell->ROmega = Cell->B[0][0]*Cell->A[0][0] +
                   Cell->B[1][0]*Cell->A[1][0] + 
                   Cell->B[2][0]*Cell->A[2][0];

    // shortest distance between opposite faces                                                                          
    for(int i=0;i<3;i++){
        Cell->rw[i]=Cell->ROmega/sqrt(Cell->A[0][i]*Cell->A[0][i]+
                                      Cell->A[1][i]*Cell->A[1][i]+ 
                                      Cell->A[2][i]*Cell->A[2][i]);
    }
    // angles ( reciprocal )
    angles_(&ang,Cell->B,Cell->Bnorm);
    for(int i=0;i<3;i++){
        Cell->rang[i]=ang[i];
    }
    /*
    printf("cell volume : %12.20f\n",omega);
    printf("%f %f %f\n",Cell->B[0][0],Cell->B[1][0],Cell->B[2][0]);
    printf("%f %f %f\n",Cell->B[0][1],Cell->B[1][1],Cell->B[2][1]);
    printf("%f %f %f\n",Cell->B[0][2],Cell->B[1][2],Cell->B[2][2]);
    */
}

void angles_(double *ang,double basis[3][3], double norm[3]){
    ang[0] = basis[0][2] * basis[0][1] + basis[1][2] * basis[1][1] + basis[2][2] * basis[2][1];
    ang[1] = basis[0][0] * basis[0][2] + basis[1][0] * basis[1][0] + basis[2][0] * basis[2][2];
    ang[2] = basis[0][0] * basis[0][1] + basis[1][0] * basis[1][1] + basis[2][0] * basis[2][1];

    ang[0] /= ( norm[2] * norm[1] );
    ang[1] /= ( norm[0] * norm[2] );
    ang[2] /= ( norm[0] * norm[1] );

    for(int i=0;i<3;i++){
        ang[i] = acos ( ang[i] ) * radian;
    }

}


/*******************************************************************************
 \brief
 transform a set of vectors from cartesian coordinates to                     
 ) direct lattice      (BASIS must be equal to B reciprocal lattice)          
 ) reciprocal lattice  (BASIS must be equal to A direct lattice)              
 \author                                                                      
 gK (VASP)                                                                    
 \note                                                                        
 adapted from VASP                                                            
 \param[in] NMAX dimension of vectors VX , VY , VZ                            
 \param[in,out] VX , VY , VZ vectors being transformed                        
 \param[in] BASIS basis vector ( direct or reciprocal lattice )               
*******************************************************************************/

void kardir (int n, double *vx, double *vy, double *vz , double basis[3][3])
{
    double v1, v2, v3;
    for (int i=0;i<n;i++)
    {
//        printf("inside kardir %d %f",i,vx[i]);
        v1=vx[i]*basis[0][0]+vy[i]*basis[0][1]+vz[i]*basis[0][2];
        v2=vx[i]*basis[1][0]+vy[i]*basis[1][1]+vz[i]*basis[1][2];
        v3=vx[i]*basis[2][0]+vy[i]*basis[2][1]+vz[i]*basis[2][2];
        vx[i]=v1;
        vy[i]=v2;
        vz[i]=v3;
//        printf(" %f\n",vx[i]);
    }
}

void dirkar(int n, double *vx, double *vy, double *vz , double basis[3][3])
{
    double v1, v2, v3;
    for (int i=0;i<n;i++)
    {
        v1=vx[i]*basis[0][0]+vy[i]*basis[1][0]+vz[i]*basis[2][0];
        v2=vx[i]*basis[0][1]+vy[i]*basis[1][1]+vz[i]*basis[2][1];
        v3=vx[i]*basis[0][2]+vy[i]*basis[1][2]+vz[i]*basis[2][2];
        vx[i]=v1;
        vy[i]=v2;
        vz[i]=v3;
    }

}


void info_cell_(char* label     , double basis[3][3],\
                double norm[3]  , double width[3]   ,\
                double angles[3], double volume ){

    char *tvec;
    tvec= ((strcmp(label,"Direct") == 0 ) ) ? "_vector " : "*_vector";
    LSEPARATOR;
    printf("%16s basis : \n",label);
    LSEPARATOR;
    for(int j=0;j<3;j++){
        printf("%c%s %14s",j+LOWERCASE,tvec,"=");
        for (int i=0;i<3;i++){
            printf("%12.4f",basis[j][i]);
        }
        putchar('\n');
    }
    putchar('\n');
    printf("cell param.            =");
    for (int i=0;i<3;i++){
        printf("%12.4f",norm[i]);
    }
    printf("\nperpend. width         =");
    for (int i=0;i<3;i++){
        printf("%12.4f",width[i]);
    }
    printf("\nangles                 =");
    for (int i=0;i<3;i++){
        printf("%12.4f",angles[i]);
    }
    printf("\nvolume                 =%12.4f\n",volume);
}

void info_simuCell(){

    if (ionode) {
        info_cell_("Direct",simuCell.A,simuCell.Anorm,simuCell.w,simuCell.ang,simuCell.Omega);
        info_cell_("Reciprocal",simuCell.B,simuCell.Bnorm,simuCell.rw,simuCell.rang,simuCell.ROmega);
        LSEPARATOR;
    }
}




