#include <stdio.h>
#include <math.h>
#include "cell.h"
#include "tools.h"

void lattice(CELL * Cell)
{
    Cell->G[0][0] = Cell->A[0][0]*Cell->A[0][0] + Cell->A[1][0]*Cell->A[1][0]  + Cell->A[2][0]*Cell->A[2][0] ;
    Cell->G[1][0] = Cell->A[1][0]*Cell->A[1][1] + Cell->A[2][0]*Cell->A[2][1];
    Cell->G[2][0] = Cell->A[2][0]*Cell->A[2][2];
    Cell->G[1][1] = Cell->A[1][1]*Cell->A[1][1] + Cell->A[2][2]*Cell->A[2][2];
    Cell->G[2][1] = Cell->A[2][1]*Cell->A[2][2];
    Cell->G[0][2] = Cell->G[2][0];
    Cell->G[1][2] = Cell->G[2][1];
    Cell->G[2][2] = Cell->A[2][2]*Cell->A[2][2];

    double u1[3],u2[3],tmp[3];

    expro(Cell->B[0],Cell->A[1],Cell->A[2]);
    expro(Cell->B[1],Cell->A[2],Cell->A[0]);
    expro(Cell->B[2],Cell->A[0],Cell->A[1]);

    
    double omega;
    omega = Cell->B[0][0]*Cell->A[0][0]+Cell->B[1][0]*Cell->A[1][0] + Cell->B[2][0]*Cell->A[2][0];
    Cell->Omega=omega;

    for(int i=0;i<3;i++) {
        for(int j=0;j<3;j++) {
            Cell->B[i][j]=Cell->B[i][j]/omega ;
        }
    }

    for (int i=0;i<3;i++) {
        for(int j=0;j<3;j++) {
            Cell->Anorm[i]+=Cell->A[j][i]*Cell->A[j][i];
            Cell->Bnorm[i]+=Cell->B[j][i]*Cell->B[j][i];
        }
        Cell->Anorm[i]=sqrt(Cell->Anorm[i]);
        Cell->Bnorm[i]=sqrt(Cell->Bnorm[i]);
    }


    printf("%f %f %f\n",Cell->B[0][0],Cell->B[1][0],Cell->B[2][0]);
    printf("%f %f %f\n",Cell->B[0][1],Cell->B[1][1],Cell->B[2][1]);
    printf("%f %f %f\n",Cell->B[0][2],Cell->B[1][2],Cell->B[2][2]);
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
