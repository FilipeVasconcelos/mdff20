#ifndef CELL_H
#define CELL_H

struct CELL simu_cell,*psimu_cell;
typedef struct CELL {
    double A[3][3];               // direct basis vector 
    double B[3][3];               // reciprocal basis vectors 
    double G[3][3];               // metric tensor (A^T A)
    double Anorm[3];              // norm of direct basis vectors
    double Bnorm[3];              // norm of reciprocal basis vectors 
    double Omega;                 // volume ( direct )
    double ROmega;                // volume ( reciprocal )
    double wa, wb,wc;             // perpendicular width (direct)
    double alpha, beta, gamma;    // angles ( direct ) 
    double rwa, rwb, rwc;         // perpendicular width (reciprocal)
    double ralpha, rbeta, rgamma; // angles ( reciprocal )    
} CELL;


void lattice(CELL * Cell);
void kardir (int n, double *vx, double *vy, double *vz , double basis[3][3]);
void dirkar (int n, double *vx, double *vy, double *vz , double basis[3][3]);
#endif /* CELL_H */
