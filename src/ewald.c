#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "field.h"
#include "tensors.h"
#include "global.h"
#include "cell.h"
#include "verlet.h"
#include "ewald.h"
#include "pbc.h"

/******************************************************************************/
void get_monopoles(){
}

/******************************************************************************/
void get_dipoles(){
}

/******************************************************************************/
void get_quadrupoles(){
}


/******************************************************************************/
void alloc_multipole(){
}


/******************************************************************************/
/* Ewald Summation */
void multipole_ES(double *q, double (*mu)[3], double (*theta)[3][3]){

    double u_dir, u_rec;
    double (*ef_dir)[3], (*ef_rec)[3];
    double (*efg_dir)[3][3], (*efg_rec)[3][3];
    double *fx_dir, *fy_dir, *fz_dir;
    double *fx_rec, *fy_rec, *fz_rec;
    double tau_dir[3][3],tau_rec[3][3];

    /*************************************************************/
    /*                DIRECT SPACE                               */
    /*************************************************************/
    ef_dir=malloc(nion*sizeof(*ef_dir));
    efg_dir=malloc(nion*sizeof(*efg_dir));
    fx_dir=malloc(nion*sizeof(*fx_dir));
    fy_dir=malloc(nion*sizeof(*fy_dir));
    fz_dir=malloc(nion*sizeof(*fz_dir));
    
    multipole_ES_dir(q, mu, theta, u_dir, ef_dir, efg_dir, fx_dir, fy_dir, fz_dir, tau_dir);

    /*************************************************************/
    /*                RECIPROCAL SPACE                           */
    /*************************************************************/
    ef_rec=malloc(nion*sizeof(*ef_rec));
    efg_rec=malloc(nion*sizeof(*efg_rec));
    fx_rec=malloc(nion*sizeof(*fx_rec));
    fy_rec=malloc(nion*sizeof(*fy_rec));
    fz_rec=malloc(nion*sizeof(*fz_rec));

    multipole_ES_rec(q, mu, theta, u_rec, ef_rec, efg_rec, fx_rec, fy_rec, fz_rec, tau_rec);

    free(ef_dir);
    free(ef_rec);
    free(efg_dir);
    free(efg_rec);
    free(fx_dir);free(fy_dir);free(fz_dir);
    free(fx_rec);free(fy_rec);free(fz_rec);
}


/******************************************************************************/
/* Multipole expansion of the Ewald sum  in direct space                      */
void multipole_ES_dir(double *q, double (*mu)[3], double (*theta)[3][3], 
                                double u_dir, double (*ef_dir)[3], double (*efg_dir)[3][3], 
                                double *fx_dir, double *fy_dir, double *fz_dir , double tau_dir[3][3]){

    double rxi,ryi,rzi;
    double rxij,ryij,rzij;
    double rij[3];
    double fij[3];
    int ja,jb,je;
    TENSOR_RK0 T0;
    TENSOR_RK1 T1;
    TENSOR_RK2 T2;
    TENSOR_RK3 T3;
    TENSOR_RK4 T4;
    double qi;
    double mui[3];
    double thetai[3][3];
    double qj;
    double muj[3];
    double thetaj[3][3];

    double expon;
    double qij;
    double d , d2 , d3  , d5 , d7 , d9;
    double dm1 , dm3 , dm5  , dm7 , dm9 , dm11;
    double F0, F1, F2, F3, F4, F5;
    double F1_dm3 , F1d_dm3 , F1d2_dm3 , F2_dm5 , F2d_dm5 , F2d2_dm5 , F3_dm7 , F4_dm9 , F5_dm11;

    double cutsq;
    double alpha2, alpha3, alpha5, alpha7, alpha9;


    //  few constants                                         
    cutsq  = verlet_coul->cut * verlet_coul->cut; //cutlongrange
    alpha2 = alphaES * alphaES;
    alpha3 = alpha2  * alphaES;
    alpha5 = alpha3  * alpha2;
    alpha7 = alpha5  * alpha2;
    alpha9 = alpha7  * alpha2;

    /*************************************** 
            cartesian to direct                    
     ***************************************/
    kardir ( nion , rx , ry , rz , simuCell.B ) ;

    for(int ia=atomDec.iaStart;ia<atomDec.iaEnd;ia++) {
        rxi = rx[ia];
        ryi = ry[ia];
        rzi = rz[ia];
        qi  = q[ia];
        for(int i=0;i<3;i++){
            mui[i] = mu[ia][i];
            for(int j=0;j<3;j++){
                thetai[i][j] = theta[ia][i][j];
            }
        }
        if (lverletL) {
            jb=verlet_coul->point[ia];
            je=verlet_coul->point[ia+1];
        }
        else{
            jb = 0;
            je = nion;
        }
    
        for (int j1=jb;j1<je;j1++) {
            if (lverletL){ 
                ja=verlet_nb->list[j1];
            }
            else{
                ja=j1;
            }
            if ( ! ( (( ja > ia ) && !lverletL ) || (( ja !=ia ) && lverletL) ) ) continue;

            qj  = q[ja];
            for(int i=0;i<3;i++){
                muj[i] = mu[ja][i];
                for(int j=0;j<3;j++){
                    thetaj[i][j] = theta[ja][i][j];
                }
            }
            for (int i=0;i<3;i++){
                fij[i] = 0.0;
            }
            qij = qi * qj ;
            rij[0]= rxi - rx[ja];
            rij[1]= ryi - ry[ja];
            rij[2]= rzi - rz[ja];
            pbc(&rij[0],&rij[1],&rij[2]);
            d2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
            if ( d2 > cutsq ) continue ;
            d    = sqrt(d2);
            d3   = d2 * d  ;
            d5   = d3 * d2 ;
            d7   = d5 * d2 ;
            d9   = d7 * d2 ;
            dm1  = 1.0 / d ; 
            dm3  = dm1 / d2;
            dm5  = dm3 / d2;
            dm7  = dm5 / d2;
            dm9  = dm7 / d2;
            dm11 = dm9 / d2;


            expon = exp( - alpha2 * d2 )/ piroot;              
            F0    = erfc( alphaES * d );                           
            F1    = F0 +  2.0 * alphaES * d  * expon;            
            F2    = F1 +  4.0 * alpha3  * d3 * expon / 3.0;   
            F3    = F2 +  8.0 * alpha5  * d5 * expon / 15.0;
            F4    = F3 + 16.0 * alpha7  * d7 * expon / 105.0; 
            F5    = F4 + 32.0 * alpha9  * d9 * expon / 945.0; 


            /*****************************************/
            /* multipole interaction tensor rank = 0 */
            /*****************************************/
            T0.sca = dm1 * F0;

            /*****************************************/
            /* multipole interaction tensor rank = 1 */
            /*****************************************/
            for (int i=0; i<3;i++){
                T1.a[i] = - rij[i] * dm3 * F1;
            }

            /*****************************************/
  	    /*  multipole interaction tensor rank = 2*/                             
            /*  nb of components = 9 => reduced = 6  */                              
            /*  + damping                            */                             
            /*****************************************/
            F2_dm5 = 3.0 * F2 * dm5;                                        
            F1_dm3 = F1  * dm3;                                                 
            for(int j=0; j<3;j++){
		 for(int k=0; k<3;k++){
                    T2.ab[j][k]= rij[j] * rij[k] * F2_dm5;
                    if (j == k ) T2.ab[j][k] += -F1_dm3; 
		}
	    }
            /*******************************************/
 	    /*   multipole interaction tensor rank = 3 */
            /*   nb of components = 27 => reduced = 10 */
            /*******************************************/
	    F3_dm7 = F3 * dm7 * 15.0;
 	    //T3.abc = 0.0;
	    for(int i=0; i<3;i++){
	    	for(int j=0; j<3;j++){
	    		for(int k=0; k<3;k++){
                            T3.abc[i][j][k]=- rij[i] * rij[j] * rij[k] * F3_dm7;
                            if ( i == j ) T3.abc[i][j][k]+= rij[k] * F2_dm5;
                            if ( i == k ) T3.abc[i][j][k]+= rij[j] * F2_dm5;
                            if ( j == k ) T3.abc[i][j][k]+= rij[i] * F2_dm5;
			}
		}
	    }

            /*******************************************/
            /*   multipole interaction tensor rank = 4 */                                                
            /*   nb of components = 81 => reduced = 15 */                                                
            /*******************************************/
            F4_dm9 = dm9 * F4 * 105.0;
	    for(int i=0; i<3;i++){
                for(int j=0; j<3;j++){
                    for(int k=0; k<3;k++){
                        for(int l=0; l<3;l++){
                            T4.abcd[i][j][k][l] = rij[i] * rij[j] * rij[k] * rij[l] * F4_dm9;
                            if ( k == l ) T4.abcd[i][j][k][l] += -rij[i]*rij[j] * F3_dm7;
                            if ( j == l ) T4.abcd[i][j][k][l] += -rij[i]*rij[k] * F3_dm7;
                            if ( j == k ) T4.abcd[i][j][k][l] += -rij[i]*rij[l] * F3_dm7;
                            if ( i == l ) T4.abcd[i][j][k][l] += -rij[j]*rij[k] * F3_dm7;
                            if ( i == k ) T4.abcd[i][j][k][l] += -rij[j]*rij[l] * F3_dm7;
                            if ( i == j ) T4.abcd[i][j][k][l] += -rij[k]*rij[l] * F3_dm7;
                            if ( ( i == j ) && ( k == l ) ) T4.abcd[i][j][k][l] += F2_dm5;
                            if ( ( i == k ) && ( j == l ) ) T4.abcd[i][j][k][l] += F2_dm5;
                            if ( ( i == l ) && ( j == k ) ) T4.abcd[i][j][k][l] += F2_dm5;
                        }
                    }
                }
            }
	
            /*
        	note :                                                                    
        	les termes faisant intervenir les tenseurs d'interactions : T0,T2,T4...   
	        sont symétriques par changement de direction de l'interaction rij => rji  
                alors que T1,T3,T5 ... ne le sont pas.                                    
                en pratique : les sommes sur i ou j change de signe pour T1,T3,T5         
	                                                                           
            */ 
            /*            
               ===========================================================               
	                      charge-charge interaction                                
               ===========================================================  
            */
	    
            /*  energy */
            u_dir += qij * T0.sca;

            /* electric field */
            for(int i=0; i<3;i++){
	        ef_dir [ia][i] += -qj * T1.a[i] ;
                ef_dir [ja][i] +=  qi * T1.a[i] ;
            }
            /* electric field gradient */
            for(int i=0; i<3;i++){
                for(int j=0; j<3;j++){
                    efg_dir [ia][i][j] += -qj * T2.ab[i][j];
                    efg_dir [ja][i][j] += -qi * T2.ab[i][j];
                }
            }
            /* forces */
            for(int i=0; i<3;i++){
                fij[i] += qij * T1.a[i];
            }




            /* TOTAL FORCES */
            fx_dir[ia] = fx_dir[ia] - fij[0];
            fy_dir[ia] = fy_dir[ia] - fij[1];
            fz_dir[ia] = fz_dir[ia] - fij[2];
            fx_dir[ja] = fx_dir[ja] + fij[0];
            fy_dir[ja] = fy_dir[ja] + fij[1];
            fz_dir[ja] = fz_dir[ja] + fij[2];
            
            /* TOTAL STRESS TENSOR */
            for(int i=0; i<3;i++){    
                for(int j=0; j<3;j++){
                    tau_dir[i][j] += - (rij[i] * fij[j] + rij[j] * fij[i]);
                }
            }

        } /* ja */

    } /* ia */



    /*************************************** 
            direct to cartesian                   
     ***************************************/
    dirkar ( nion , rx , ry , rz , simuCell.A ) ;
}

/******************************************************************************/
/* Multipole expansion of the Ewald sum  in reciprocal space                      */
void multipole_ES_rec(double *q, double (*mu)[3], double (*theta)[3][3],
                      double u_rec  , double (*ef_rec)[3], double (*efg_rec)[3][3],
                      double *fx_rec, double *fy_rec, double *fz_rec , double tau_rec[3][3]){
}

