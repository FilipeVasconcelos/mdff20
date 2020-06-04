#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "config.h"
#include "field.h"
#include "tensors.h"
#include "global.h"
#include "cell.h"
#include "verlet.h"
#include "ewald.h"
#include "pbc.h"
#include "tools.h"
#include "math_mdff.h"
#include "kspace.h"
#include "functions.h"



void set_autoES(){

    double rcut;
    double eps;
    double tol,alpha,tol1;

    rcut =  0.5 * dmin_arr(simuCell.w,3); 
    printf("rcut %f min %f  \n",rcut,dmin_arr(simuCell.w,3));
    if ( cutlongrange < rcut ) {
        printf("WARNING : cutlongrange will be changed according to simuCell.W[]\n");
        cutlongrange=rcut;
    }
    eps=dmin(fabs(epsw),0.5);
    tol=sqrt(fabs(log(eps*cutlongrange)));
    alpha=sqrt(abs(log(eps*cutlongrange*tol)))/cutlongrange;
    tol1=sqrt(-log(eps*cutlongrange*(2.0*tol*alpha)*(2.0*tol*alpha)));

    printf("automatic Ewald Sum\n");    
    alphaES = alpha;
    for(int k=0;k<3;k++){
        kES[k]=(int) nint(0.25+simuCell.Anorm[k]*alpha*tol1/PI);
    }

//    printf("eps %e tol %e  tol1 %e %d\n",eps,tol,tol1,((int) nint(0.25+simuCell.Anorm[0]*alpha*tol1/PI)));
}

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
void multipole_ES(double *q, double (*mu)[3], double (*theta)[3][3],double *u){

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
    
    multipole_ES_dir(q, mu, theta, &u_dir, ef_dir, efg_dir, fx_dir, fy_dir, fz_dir, tau_dir);

    /*************************************************************/
    /*                RECIPROCAL SPACE                           */
    /*************************************************************/
    ef_rec=malloc(nion*sizeof(*ef_rec));
    efg_rec=malloc(nion*sizeof(*efg_rec));
    fx_rec=malloc(nion*sizeof(*fx_rec));
    fy_rec=malloc(nion*sizeof(*fy_rec));
    fz_rec=malloc(nion*sizeof(*fz_rec));

    multipole_ES_rec(q, mu, theta, &u_rec, ef_rec, efg_rec, fx_rec, fy_rec, fz_rec, tau_rec);


    double tpi_V,tpi_3V;
    double qsq,selfa;
    double u_self_qq;
    double u_surf_qq;
    double qt[3];

    qsq=0.0;
    qt[0]=0;    qt[1]=0;    qt[2]=0;
    for(int ia=0;ia<nion;ia++){
        qt[0]+=q[ia]*rx[ia];
        qt[1]+=q[ia]*ry[ia];
        qt[2]+=q[ia]*rz[ia];
        qsq += q[ia];
    }
    tpi_V  = TPI    * simuCell.inveOmega;  // 2pi / V 
    tpi_3V = tpi_V  / 3.0;                 // 2pi / 3V

    u_surf_qq = qt[0]*qt[0] + qt[1]*qt[1] + qt[2]*qt[2]; 
    u_surf_qq *= tpi_3V; 

    selfa  = alphaES / piroot ;
    u_self_qq   =  - selfa  * qsq;

    //colombic potential energy
    printf("u_dir %e\n",u_dir);
    printf("u_rec %e\n",u_rec);
    printf("u_self_qq %e\n",u_self_qq);
    printf("u_surf_qq %e\n",u_surf_qq);
    *u= (u_dir + u_rec + u_self_qq + u_surf_qq )*coul_unit;


    for(int ia=0;ia<nion;ia++){
        fx[ia] = fx_dir[ia] + fx_rec[ia];
        fy[ia] = fy_dir[ia] + fy_rec[ia];
        fz[ia] = fz_dir[ia] + fz_rec[ia];
        printf("ia %d fx_dir %15.8e\n",ia,fx_dir[ia]);
        printf("ia %d fx_rec %15.8e\n",ia,fx_rec[ia]);
        printf("ia %d fx %15.8e\n",ia,fx[ia]);
    }

    sample_config(0);

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
                      double *u_dir, double (*ef_dir)[3], double (*efg_dir)[3][3], 
                      double *fx_dir, double *fy_dir, double *fz_dir , double tau_dir[3][3]){

    printf("inside multipole_ES_dir\n");
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
    double uu = 0;

    //  few constants                                         
    if (lverletL) {
        cutsq  = verlet_coul->cut * verlet_coul->cut; //cutlongrange
    }
    else{
        cutsq = cutlongrange*cutlongrange;
    }

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
        printf("ia %d\n",ia);
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
            if ( ( (( ja <= ia ) && !lverletL ) || (( ja == ia ) && lverletL) ) ) continue;

            qj  = q[ja];
            for(int i=0;i<3;i++){
                muj[i] = mu[ja][i];
                fij[i] = 0.0;
                for(int j=0;j<3;j++){
                    thetaj[i][j] = theta[ja][i][j];
                }
            }
            qij = qi * qj ;
            rij[0]= rxi - rx[ja];
            rij[1]= ryi - ry[ja];
            rij[2]= rzi - rz[ja];
            pbc(&rij[0],&rij[1],&rij[2]);
            d2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
            printf("test continue %f %f\n",d2,cutsq);
            if ( d2 > cutsq ) continue ;
            printf("passé !\n");
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
            //F0    = erfc( alphaES * d );
            F0    = errfc( alphaES * d );
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
	    
            /*  potential energy */
            uu += qij * T0.sca;
            printf("here !!! %f %f %f %f %f\n",d2,cutsq,qij,T0.sca,uu);

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
            fx_dir[ia] += -fij[0];
            fy_dir[ia] += -fij[1];
            fz_dir[ia] += -fij[2];
            fx_dir[ja] +=  fij[0];
            fy_dir[ja] +=  fij[1];
            fz_dir[ja] +=  fij[2];
            
            /* TOTAL STRESS TENSOR */
            for(int i=0; i<3;i++){    
                for(int j=0; j<3;j++){
                    tau_dir[i][j] += - (rij[i] * fij[j] + rij[j] * fij[i]);
                }
            }

        } /* ja */

    } /* ia */

    *u_dir = uu; 

    /*************************************** 
            direct to cartesian                   
     ***************************************/
    dirkar ( nion , rx , ry , rz , simuCell.A ) ;
}

/******************************************************************************/
/* Multipole expansion of the Ewald sum  in reciprocal space                  */
/******************************************************************************/
void multipole_ES_rec(double *q, double (*mu)[3], double (*theta)[3][3],
                      double *u_rec  , double (*ef_rec)[3], double (*efg_rec)[3][3],
                      double *fx_rec, double *fy_rec, double *fz_rec , double tau_rec[3][3]){

    double uu;
    double kx,ky,kz,kk,Ak,kcoe;
    double rxi,ryi,rzi;
    double qi,k_dot_r;
    double *ckr,*skr;
    double rhonk_R,rhonk_I;
    double recarg;

    ckr=malloc(nion*sizeof(*ckr));
    skr=malloc(nion*sizeof(*skr));

    uu=0;

    for (int ik=kcoul.kptDec.iaStart;ik<kcoul.kptDec.iaEnd;ik++){
        if (kcoul.kk[ik] == 0.0) continue;
        kx   = kcoul.kx[ik];
        ky   = kcoul.ky[ik];
        kz   = kcoul.kz[ik];
        kk   = kcoul.kk[ik];
        Ak   = kcoul.Ak[ik];
        kcoe = kcoul.kcoe[ik];

        rhonk_R = 0.0;
        rhonk_I = 0.0;

        for (int ia=0;ia<nion;ia++){
            rxi = rx[ia];
            ryi = ry[ia];
            rzi = rz[ia];
            qi  = q[ia];
            k_dot_r  = ( kx * rxi + ky * ryi + kz * rzi );
            ckr[ia]  = cos(k_dot_r);
            skr[ia]  = sin(k_dot_r);
            rhonk_R  += qi * ckr[ia];
            rhonk_I  += qi * skr[ia];
//            printf("%d %d %e %e %e %e\n",ik,ia,rhonk_R,qi * ckria,rhonk_I,qi * skria);
            
        } /* sum to get charge density */
        /* potential energy */                                       
        uu  += (rhonk_R*rhonk_R + rhonk_I*rhonk_I) * Ak ;                                  

        /* second sum */
        for (int ia=0;ia<nion;ia++){
            qi=q[ia];
            recarg  = Ak * (rhonk_I*ckr[ia] - rhonk_R*skr[ia]);
            fx_rec[ia] += -qi * kx * recarg;
            fy_rec[ia] += -qi * ky * recarg;
            fz_rec[ia] += -qi * kz * recarg;
        }
        //printf("%d %e %e %e %e\n",ik,Ak,rhonk_R,rhonk_I,uu);

    }

    double tpiV = TPI * simuCell.inveOmega;
    double fpiV = 2.0 * tpiV;
    /* half mesh */
    uu = 2.0 * uu;
    for (int ia=0;ia<nion;ia++){
        fx_rec[ia] *= 2.0 * fpiV;
        fy_rec[ia] *= 2.0 * fpiV;
        fz_rec[ia] *= 2.0 * fpiV;
    }

    /* remark on unit
    1/(4*pi*epislon_0) = 1 => epsilon_0 = 1/4pi
    */
    *u_rec =   uu * tpiV;


    free(ckr);free(skr);

}

