#ifndef CELL_H
#define CELL_H

struct CELL simuCell;
typedef struct CELL {
    double A[3][3];               // direct basis vector 
    double B[3][3];               // reciprocal basis vectors 
    double G[3][3];               // metric tensor (A^T A)
    double Anorm[3];              // norm of direct basis vectors
    double Bnorm[3];              // norm of reciprocal basis vectors 
    double Omega;                 // Ω : volume ( direct )
    double inveOmega;             // 1/Ω (direct)
    double inveOmegaU;            // 1/Ω (direct) in press_unit
    double ROmega;                // volume ( reciprocal )
    double w[3];                  // perpendicular width (direct)
    double ang[3];                // angles ( direct ) 
    double rw[3];                 // perpendicular width (reciprocal)
    double rang[3];               // angles ( reciprocal )    
} CELL;


void lattice(CELL * Cell);
void kardir (int n, double *vx, double *vy, double *vz , double basis[3][3]);
void dirkar (int n, double *vx, double *vy, double *vz , double basis[3][3]);
void info_simuCell();
void angles_();
#endif /* CELL_H */
