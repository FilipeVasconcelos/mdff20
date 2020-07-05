#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "coulombic.h"
#include "global.h"
#include "config.h"
#include "field.h"
#include "tensors.h"
#include "global.h"
#include "cell.h"
#include "verlet.h"
#include "ewald.h"
#include "pbc.h"
#include "timing.h"
#include "tools.h"
#include "math_mdff.h"
#include "kspace.h"
#include "functions.h"
#include "thermo.h"
#include "tt_damp.h"

//#define DEBUG_EWALD
//#define DEBUG_EWALD_DIR

/******************************************************************************/
void set_autoES(){

    double rcut;
    double eps;
    double tol,alpha,tol1;

    rcut =  0.5 * dmin_arr(simuCell.w,3); 
//    printf("rcut %f min %f  \n",rcut,dmin_arr(simuCell.w,3));
    if ( cutlongrange < rcut ) {
        printf("WARNING : cutlongrange will be changed according to simuCell.W[]\n");
        cutlongrange=rcut;
    }
    eps=dmin(fabs(epsw),0.5);
    tol=sqrt(fabs(log(eps*cutlongrange)));
    alpha=sqrt(abs(log(eps*cutlongrange*tol)))/cutlongrange;
    tol1=sqrt(-log(eps*cutlongrange*(2.0*tol*alpha)*(2.0*tol*alpha)));

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
void get_quadrupoles(){
}


/******************************************************************************/
void alloc_multipole(){
}


/******************************************************************************/
/* Ewald Summation */
void multipole_ES(double *q, double (*mu)[3], double (*theta)[3][3],double *u, double *pvir, double tau[3][3],
                             double (*ef)[3], double (*efg)[3][3], bool lq, bool ld, bool lcouldamp){

    double u_dir, u_rec;
    double (*ef_dir)[3], (*ef_rec)[3], (*ef_self)[3];
    double (*efg_dir)[3][3], (*efg_rec)[3][3], (*efg_self)[3][3];
    double *fx_dir, *fy_dir, *fz_dir;
    double *fx_rec, *fy_rec, *fz_rec;
    double tau_dir[3][3],tau_rec[3][3];

    statime(18);

    u_dir=0.0;
    u_rec=0.0;
    /*************************************************************/
    /*                DIRECT SPACE                               */
    /*************************************************************/
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            tau_dir[i][j]=0.0;
        }
    }
    ef_dir=malloc(nion*sizeof(*ef_dir));
    efg_dir=malloc(nion*sizeof(*efg_dir));
    fx_dir=malloc(nion*sizeof(*fx_dir));
    fy_dir=malloc(nion*sizeof(*fy_dir));
    fz_dir=malloc(nion*sizeof(*fz_dir));
    for (int ia=0;ia<nion;ia++){
        fx_dir[ia]=0.0;
        fy_dir[ia]=0.0;
        fz_dir[ia]=0.0;
        for (int i=0;i<3;i++){
            ef_dir[ia][i]=0.0;
            for (int j=0;j<3;j++){
                efg_dir[ia][i][j]=0.0;
            }
        }
    }
/*
#ifdef DEBUG_EWALD
    for (int i=0;i<3;i++){
        printf("ef_dir  %f\n",ef_dir[0][i]);
        printf("ef_dir  %f\n",ef_dir[1][i]);
        printf("\n");
    }
#endif 
*/
    multipole_ES_dir(q, mu, theta, &u_dir, ef_dir, efg_dir, fx_dir, fy_dir, fz_dir, tau_dir, lq, ld ,lcouldamp);
/*
#ifdef DEBUG_EWALD
    for (int i=0;i<3;i++){
        printf("ef_dir  %f\n",ef_dir[0][i]);
        printf("ef_dir  %f\n",ef_dir[1][i]);
        printf("\n");
    }
#endif 
*/
    statime(19);
    mestime(&ewaldDirCPUtime,19,18);

    /*************************************************************/
    /*                RECIPROCAL SPACE                           */
    /*************************************************************/
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            tau_rec[i][j]=0.0;
        }
    }
    ef_rec=malloc(nion*sizeof(*ef_rec));
    efg_rec=malloc(nion*sizeof(*efg_rec));
    fx_rec=malloc(nion*sizeof(*fx_rec));
    fy_rec=malloc(nion*sizeof(*fy_rec));
    fz_rec=malloc(nion*sizeof(*fz_rec));
    for (int ia=0;ia<nion;ia++){
        fx_rec[ia]=0.0;
        fy_rec[ia]=0.0;
        fz_rec[ia]=0.0;
        for (int i=0;i<3;i++){
            ef_rec[ia][i]=0.0;
            for (int j=0;j<3;j++){
                efg_rec[ia][i][j]=0.0;
            }
        }
    }

    multipole_ES_rec(q, mu, theta, &u_rec, ef_rec, efg_rec, fx_rec, fy_rec, fz_rec, tau_rec, lq, ld);
    
    statime(20);
    mestime(&ewaldRecCPUtime,20,19);

    /*************************************************************/
    /*                SELF AND SURFACE                           */
    /*************************************************************/
    double alpha2;
    double qsq,musq,selfa,selfa2;
    double u_self_qq,u_self_dd,u_self;
    double qt[3];

    qsq=0.0;
    qt[0]=0.0;    qt[1]=0.0;    qt[2]=0.0;
    musq=0.0;
    for(int ia=0;ia<nion;ia++){
        qt[0]+=q[ia]*rx[ia];
        qt[1]+=q[ia]*ry[ia];
        qt[2]+=q[ia]*rz[ia];
        qsq  += q[ia]*q[ia];
        musq += mu[ia][0]*mu[ia][0] + mu[ia][1]*mu[ia][1] + mu[ia][2]*mu[ia][2];
    }

    alpha2 = alphaES * alphaES;
    selfa  = alphaES / piroot ;
    selfa2 = 2.0 * selfa * alpha2 / 3.0;
    u_self_qq = -selfa * qsq;
    u_self_dd = -selfa2 * musq;
    u_self = u_self_qq + u_self_dd;
    
    ef_self=malloc(nion*sizeof(*ef_self));
    efg_self=malloc(nion*sizeof(*efg_self));
    for (int ia=0;ia<nion;ia++){
        for (int i=0;i<3;i++){
            ef_self[ia][i]= 2.0 * selfa2 *mu[ia][i];
            efg_self[ia][i][i]=- 2.0 * selfa2 * q[ia];
            for (int j=0;j<3;j++){
                if (j!=i) efg_self[ia][i][j]=0.0;
            }
        }
    }

    /* ------------------------------------------ */
    /*  total quantities and electrostatic units  */
    /* ------------------------------------------ */
    /* electrostatic potential energy */    
    *u= (u_dir + u_rec + u_self ) * coul_unit;

    for (int ia=0;ia<nion;ia++){
        /* --------------------------*/
        /* electric field at ions    */
        /* --------------------------*/
        for (int i=0;i<3;i++){
            ef[ia][i]=(ef_dir[ia][i]+ef_rec[ia][i]+ef_self[ia][i]);
#ifdef DEBUG_EWALD
            printf("EF : %d %d dir: %e rec: %e self : %e tot : %e\n",ia,i,ef_dir[ia][i],ef_rec[ia][i],ef_self[ia][i],ef[ia][i]);
#endif 
        /* ------------------------- */
        /* electric field gradient   */
        /* ------------------------- */
            for (int j=0;j<3;j++){
                efg[ia][i][j]=(efg_dir[ia][i][j]+efg_rec[ia][i][j]+efg_self[ia][i][j]);
#ifdef DEBUG_EWALD
        printf("EFG : %d %d %d dir: %e rec: %e self : %e tot : %e\n",ia,i,j,efg_dir[ia][i][j],efg_rec[ia][i][j],efg_self[ia][i][j],efg[ia][i][j]);
#endif 
            }
        }
        /* ---------------------------- */
        /*           forces             */
        /* ---------------------------- */
        fx[ia] += (fx_dir[ia] + fx_rec[ia]) * coul_unit;
        fy[ia] += (fy_dir[ia] + fy_rec[ia]) * coul_unit;
        fz[ia] += (fz_dir[ia] + fz_rec[ia]) * coul_unit;
    }

    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            tau[i][j] = ( tau_dir[i][j] + tau_rec[i][j] ) * coul_unit;
        }
    }
    
    *pvir=0.0;
    for (int i=0;i<3;i++){
        *pvir+= tau[i][i] * onethird ;
    }

#ifdef DEBUG_EWALD
    putchar('\n');
    for (int ia=0;ia<nion;ia++){
        printf("mu  %e %e %e\n",mu[ia][0],mu[ia][1],mu[ia][2]);
    }
    printf("u_dir  %e\n",u_dir*coul_unit);
    printf("u_rec  %e\n",u_rec*coul_unit);
    printf("u_self %e\n",u_self*coul_unit);
    printf("u_coul %e\n",*u);
    printf("fx_dir  %e %e %e \n",fx_dir[0],fy_dir[0],fz_dir[0]);
    printf("fx_rec  %e %e %e \n",fx_rec[0],fx_rec[0],fx_rec[0]);
    printf("ef_dir  %e %e %e \n",ef_dir[0][0],ef_dir[0][1],ef_dir[0][2]);
    printf("ef_rec  %e %e %e\n",ef_rec[0][0],ef_rec[0][1],ef_rec[0][2]);
    putchar('\n');
#endif
    
    free(ef_dir);
    free(ef_rec);
    free(ef_self);
    free(efg_dir);
    free(efg_rec);
    free(efg_self);
    free(fx_dir);
    free(fy_dir);
    free(fz_dir);
    free(fx_rec);
    free(fy_rec);
    free(fz_rec);
}


/******************************************************************************/
/* Multipole expansion of the Ewald sum  in direct space                      */
void multipole_ES_dir(double *q, double (*mu)[3], double (*theta)[3][3], 
                      double *u_dir, double (*ef_dir)[3], double (*efg_dir)[3][3], 
                      double *fx_dir, double *fy_dir, double *fz_dir , double tau_dir[3][3],
                      bool lqchtask, bool ldiptask, bool lcouldamp){

    double rxi,ryi,rzi;
    double rij[3];
    double fxij,fyij,fzij;
    TENSOR_RK0 T0;
    TENSOR_RK1 T1;
    TENSOR_RK2 T2;
    TENSOR_RK3 T3;
    TENSOR_RK4 T4;
    double qi;
    double mui[3];
    //double thetai[3][3];
    double qj;
    double muj[3];
    bool ldamp;
    bool qch_i,qch_j,qch_iETj,qch_iOUj;
    bool dip_i,dip_j,dip_iETj,dip_iOUj;
    //double thetaj[3][3];
    /*
    double *ef_x,*ef_y,*ef_z;
    double ef_xi,ef_yi,ef_zi;
    double ef_xj,ef_yj,ef_zj;
    ef_x=malloc(nion*sizeof(ef_x));
    ef_y=malloc(nion*sizeof(ef_y));
    ef_z=malloc(nion*sizeof(ef_z));
    double (*tefg)[3][3];
    tefg=malloc(nion*sizeof(*tefg));
    for (int ia=0;ia<nion;ia++){
        ef_x[ia]=0.0;
        ef_y[ia]=0.0;
        ef_z[ia]=0.0;
        for (int i=0;i<3;i++){
            for (int j=0;j<3;j++){
                tefg[ia][i][j]=0.0;
            }
        }
    }
    */



    double expon;
    double qij;
    double d , d2 , d3  , d5 , d7, d9;
    double dm1 , dm3 , dm5  , dm7 , dm9;
    double F0, F1, F2, F3, F4, F5;
    double F1_dm3 ,F2_dm5 , F3_dm7 , F4_dm9; 

    // damping related
    double fdamp1,fdamp2,fdampdiff1,fdampdiff2;
    double F1d1, F2d1, F1d2, F2d2;
    double F2d1_dm5,F2d2_dm5,F1d1_dm3,F1d2_dm3;

    double alpha2, alpha3, alpha5, alpha7, alpha9;
    double uu = 0;
    double uudmp = 0;
    int ia,j1,ja,jb,je,ita,jta;
    double ttau[3][3];
    for (int i=0;i<3;i++){
        mui[i]=0.0;
        muj[i]=0.0;
        for(int j=i;j<3;j++){
            ttau[i][j]=0.0;
        }
    }
    bool lqchdiptask = lqchtask && ldiptask;

    //  few constants                                         
    alpha2 = alphaES * alphaES;
    alpha3 = alpha2  * alphaES;
    alpha5 = alpha3  * alpha2;
    alpha7 = alpha5  * alpha2;
    alpha9 = alpha7  * alpha2;

    /*************************************** 
            cartesian to direct                    
     ***************************************/
    kardir ( nion , rx , ry , rz , simuCell.B ) ;

//#pragma omp parallel shared(srcutsq,rx,ry,rz,q,vx,vy,vz,fx,fy,fz,typia) \                                                                            private(ia,j1,jb,je,ja,rxi,ryi,rzi,rxij,ryij,rzij,rijsq,p1,p2,wij,fxij,fyij,fzij,ir2,d,erh,ir6,ir8,irf6,fdiff8,ir7,ir9)                                                                                                                {                                                                                                                                 #pragma omp for reduction (+:uu,ttau,fx[:nion],fy[:nion],fz[:nion],uudmp) schedule(dynamic,16)        

    #pragma omp parallel default(none) \
    shared(lrcutsq,q,mu,rx,ry,rz,typia,piroot,pol_damp_b,pol_damp_c,pol_damp_k,lpoldamping,nion,efg_dir,ef_dir,lqchdiptask,lqchtask,ldiptask,alphaES,alpha2,alpha3,alpha5,alpha7,alpha9,verlet_coul,fx_dir,fy_dir,fz_dir,ttau,uu,atomDec,lverletL,lcouldamp,uudmp) \
    private(ita,jta,ia,j1,jb,je,ja,qi,qj,qij,mui,muj,rij,rxi,ryi,rzi,d2,qch_i,dip_i,qch_j,dip_j,fxij,fyij,fzij,dip_iETj,dip_iOUj,qch_iETj,qch_iOUj,d,d3,d5,d7,d9,dm1,dm3,dm5,dm7,dm9,T0,T1,T2,T3,T4,fdamp1,fdamp2,fdampdiff1,fdampdiff2,expon,F0,F1,F2,F3,F4,F5,F1d1,F2d1, F1d2, F2d2,F1_dm3 ,F2_dm5 , F3_dm7 , F4_dm9, ldamp,F2d1_dm5,F2d2_dm5,F1d1_dm3,F1d2_dm3) 
    {
        #pragma omp for reduction (+:uu,ttau,fx_dir[:nion],fy_dir[:nion],fz_dir[:nion],ef_dir[:nion],efg_dir[:nion]) schedule(dynamic,16) 
        for(ia=atomDec.iaStart;ia<atomDec.iaEnd;ia++) {
            rxi = rx[ia];
            ryi = ry[ia];
            rzi = rz[ia];
            qi  = q[ia];
            qch_i= (qi == 0.0) ? false : true; 
            dip_i=false;
            ita=typia[ia];
            for(int i=0;i<3;i++){
                mui[i] = mu[ia][i];
                if (mui[i] != 0.0) dip_i=true;
                /*for(int j=0;j<3;j++){
                    thetai[i][j] = theta[ia][i][j];
                }*/
            }
            if (lverletL) {
                jb=verlet_coul->point[ia];
                je=verlet_coul->point[ia+1];
            }
            else{
                jb = 0;
                je = nion;
            }
        
            for (j1=jb;j1<je;j1++) {
                if (lverletL){ 
                    ja=verlet_coul->list[j1];
                }
                else{
                    ja=j1;
                }
                if ( ( (( ja <= ia ) && !lverletL ) || (( ja == ia ) && lverletL) ) ) continue;
    
                qj  = q[ja];
                qch_j= (qj == 0.0) ? false : true; 
                fxij = 0.0;
                fyij = 0.0;
                fzij = 0.0;
                dip_j=false;
                jta=typia[ja];
                for(int i=0;i<3;i++){
                    muj[i] = mu[ja][i];
                    if (muj[i] != 0.0) dip_j=true;
                    /*for(int j=0;j<3;j++){
                        thetaj[i][j] = theta[ja][i][j];
                    }*/
                }
                qij    = qi * qj ;
                rij[0] = rxi - rx[ja];
                rij[1] = ryi - ry[ja];
                rij[2] = rzi - rz[ja];
                pbc(&rij[0],&rij[1],&rij[2]);
                d2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
                if ( d2 > lrcutsq ) continue ;
                dip_iETj = dip_i && dip_j;
                dip_iOUj = dip_i || dip_j;
                qch_iETj = qch_i && qch_j;
                qch_iOUj = qch_i || qch_j;
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
   //             dm11 = dm9 / d2;
    
		// damping function 
                ldamp = ( (lpoldamping[ita][ita][jta] || lpoldamping[jta][ita][jta] ) && ( lcouldamp ) );
                if ( ldamp ) {
//                        printf("DAMPING FUNCTIONS %d %d %d %d %d %e %e %d\n",ia,ja,ita,jta,lcouldamp,pol_damp_b[ita][ita][jta],pol_damp_c[ita][ita][jta],pol_damp_k[ita][ita][jta]);
                        TT_damping_functions(pol_damp_b[ita][ita][jta],
                                             pol_damp_c[ita][ita][jta],d,&fdamp1,&fdampdiff1,
                                             pol_damp_k[ita][ita][jta]);
                        TT_damping_functions(pol_damp_b[jta][ita][jta],
                                             pol_damp_c[jta][ita][jta],d,&fdamp2,&fdampdiff2,
                                             pol_damp_k[jta][ita][jta]);             
                }
                else {
                    fdamp1     = 1.0;
                    fdamp2     = 1.0;
                    fdampdiff1 = 0.0;
                    fdampdiff2 = 0.0;
                }

                expon = exp( - alpha2 * d2 )/ piroot;              
                F0    = erfc( alphaES * d );
                //#pragma omp task
                //F0    = errfc( alphaES * d ); //buggy with omp 2/07/20
                F1    = F0 +  2.0 * alphaES * d  * expon;            
                F2    = F1 +  4.0 * alpha3  * d3 * expon / 3.0;   
                F3    = F2 +  8.0 * alpha5  * d5 * expon / 15.0;
                F4    = F3 + 16.0 * alpha7  * d7 * expon / 105.0; 
                F5    = F4 + 32.0 * alpha9  * d9 * expon / 945.0; 
   
                // damping if no damping fdamp == 1 and fdampdiff == 0                                                    
                // recursive relation (10) in J. Chem. Phys. 133, 234101 (2010)
                F1d1  = - fdamp1 + 1.0;
                F1d2  = - fdamp2 + 1.0;
                F2d1  = F1d1 + ( d / 3.0 ) * fdampdiff1;
                F2d2  = F1d2 + ( d / 3.0 ) * fdampdiff2;
                //printf("F1d1 %e F1d2 %e F2d1 %e F2d2 %e %e %e %e %e %d\n",F1d1,F1d2,F2d1,F2d2,fdamp1,fdamp2,fdampdiff1,fdampdiff2,ldamp);

                /***** set tensors to zeros */
                for(int i=0; i<3;i++){
                    for(int j=0; j<3;j++){
                        T2.ab[i][j]=0.0;
                        T2.abD1[i][j]=0.0;
                        T2.abD2[i][j]=0.0;
    		        for(int k=0; k<3;k++){
                            T3.abc[i][j][k]=0.0;
    		            for(int l=0; l<3;l++){
                                T4.abcd[i][j][k][l]=0.0;
                            }
                        }
                    }
                }
                /*****************************************/
                /* multipole interaction tensor rank = 0 */
                /*****************************************/
                T0.sca = dm1 * F0;
    
                /*****************************************/
                /* multipole interaction tensor rank = 1 */
                /*****************************************/
                for (int i=0; i<3;i++){
                    T1.a[i] = - rij[i] * dm3 * F1;
                    if ( ldamp ) {
                        T1.aD1[i] = T1.a[i] * F1d1;
                        T1.aD2[i] = T1.a[i] * F1d2;
                    }  
                }
    
                /*****************************************/
                /*  multipole interaction tensor rank = 2*/                             
                /*  nb of components = 9 => reduced = 6  */                              
                /*  + damping                            */                             
                /*****************************************/
                F2_dm5   = 3.0 * F2   * dm5;                                        
                F2d1_dm5 = 3.0 * F2d1 * dm5;
                F2d2_dm5 = 3.0 * F2d2 * dm5;

                F1_dm3   = F1   * dm3;                                                 
                F1d1_dm3 = F1d1 * dm3;                                                 
                F1d2_dm3 = F1d2 * dm3;                                                 
    

                for(int j=0; j<3;j++){
    		    for(int k=0; k<3;k++){
                        T2.ab[j][k]= rij[j] * rij[k] * F2_dm5;
                        if (j == k ) T2.ab[j][j] += -F1_dm3; 
                        if (ldamp) {
                            T2.abD1[j][k] = rij[j] * rij[k] * F2d1_dm5; 
                            T2.abD2[j][k] = rij[j] * rij[k] * F2d2_dm5;
                            if (j == k ) {
                                T2.abD1[j][j] += -F1d1_dm3; 
                                T2.abD2[j][j] += -F1d2_dm3; 
                            }

                        }
		    }
	        }

                if ( ldiptask ) {
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
                } /* if ldiptask */
	
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
                if ( (lqchtask) && (qch_iETj) ) {
#ifdef DEBUG_EWALD_DIR
                    printf("inside charge-charge interaction %d %d\n",ia,ja);
#endif
                    /*  potential energy */
                    uu += qij * T0.sca;
                    
                    // electric field 
                    for(int i=0; i<3;i++){
            	        ef_dir[ia][i] += -qj * T1.a[i];
                        ef_dir[ja][i] +=  qi * T1.a[i];
                        if ( ldamp ) {
                            ef_dir[ia][i] +=  qj * T1.aD1[i];
                            ef_dir[ia][i] += -qi * T1.aD2[i];
                        }
                    }
                    
                    // electric field gradient 
                    for(int i=0; i<3;i++){
                        for(int j=0; j<3;j++){
                            efg_dir[ia][i][j] += -qj * T2.ab[i][j];
                            efg_dir[ja][i][j] += -qi * T2.ab[i][j];
                        }
                    }

                    /* forces */
                    fxij = qij * T1.a[0];
                    fyij = qij * T1.a[1];
                    fzij = qij * T1.a[2];
                }
                /*            
                   ===========================================================               
    	                      dipole-dipole interaction                                
                   ===========================================================  
                */
                if ( ( ldiptask ) && (dip_iETj) ) {
#ifdef DEBUG_EWALD_DIR
                    printf("inside dipole-dipole interaction  %d %d\n",ia,ja);
#endif
                    /*  potential energy */
                    for(int i=0; i<3;i++){
                        for(int j=0; j<3;j++){
                            uu += - mui[i] * T2.ab[i][j] * muj[j];
                        }
                    }
                    /* electric field */
                    for(int i=0; i<3;i++){
                        for(int j=0; j<3;j++){
            	            ef_dir[ia][i] += muj[j] * T2.ab[i][j] ;
                            ef_dir[ja][i] += mui[j] * T2.ab[i][j] ;
                        }
                    }
                    /* electric field gradient */
                    for(int i=0; i<3;i++){
                        for(int j=0; j<3;j++){
                            for(int k=0; k<3;k++){
                                efg_dir[ia][i][j] +=  muj[k] * T3.abc[i][j][k];
                                efg_dir[ja][i][j] += -mui[k] * T3.abc[i][j][k];
                            }
                        }
                    }
                    /* forces */
                    for(int j=0; j<3;j++){
                        for(int k=0; k<3;k++){
                            fxij += -mui[j] * T3.abc[j][0][k] * muj[k];
                            fyij += -mui[j] * T3.abc[j][1][k] * muj[k];
                            fzij += -mui[j] * T3.abc[j][2][k] * muj[k];
                        }
                    }
                }/* ldiptask */
                /*            
                   ===========================================================               
    	                      charge-dipole interaction                                
                   ===========================================================  
                */
                if ( ( lqchdiptask ) && (dip_iOUj) ) {
#ifdef DEBUG_EWALD_DIR
                    printf("inside charge-dipole interaction %d %d\n",ia,ja);
#endif
                    /*  potential energy */
                    for(int i=0; i<3;i++){
                        uu += - qi * T1.a[i] * muj[i];
                        uu +=   qj * T1.a[i] * mui[i];
                        if (ldamp) {
                        //    printf("LDAMP1 %e %e %e %e\n",qi,T1.aD2[i],muj[i],qi *  T1.aD2[i] * muj[i]);
                        //    printf("LDAMP2 %e %e %e %e\n",qj,T1.aD1[i],mui[i],-qj *  T1.aD1[i] * mui[i]);
                            uu +=  qi *  T1.aD2[i] * muj[i]; 
                            uu += -qj *  T1.aD1[i] * mui[i]; 
                            uudmp +=  qi *  T1.aD2[i] * muj[i]; 
                            uudmp += -qj *  T1.aD1[i] * mui[i]; 
                        }
                    }
                    /* forces  
                    ----------------------------------------------------------------
                    remarque 1 : 
                    ----------------------------------------------------------------
                    pour garder la construction :                      
                        f(ia) = f(ia) - fij                                   
                        f(ja) = f(ja) + fij                                   
                    en fin de boucle uniquement,
                    nous avons ici un changement de signe sur la somme sur j malgré 
                    le fait que le tenseur d'interaction soit paire (e.g T2 )       
                    ----------------------------------------------------------------
                    remarque 2 : 
                    ----------------------------------------------------------------
                    thole functions only apply to fields not forces    
                    ----------------------------------------------------------------
                    */
                    for(int j=0; j<3;j++){
                        fxij += -qi * T2.ab[0][j] * muj[j];
                        fyij += -qi * T2.ab[1][j] * muj[j];
                        fzij += -qi * T2.ab[2][j] * muj[j];
                        fxij +=  qj * T2.ab[0][j] * mui[j];
                        fyij +=  qj * T2.ab[1][j] * mui[j];
                        fzij +=  qj * T2.ab[2][j] * mui[j];
                    }
                    if (ldamp) {
                        for(int j=0; j<3;j++){
                            fxij += -qj * T2.abD1[j][0] * mui[j];
                            fyij += -qj * T2.abD1[j][1] * mui[j];
                            fzij += -qj * T2.abD1[j][2] * mui[j];
                            fxij +=  qi * T2.abD2[j][0] * muj[j];
                            fyij +=  qi * T2.abD2[j][1] * muj[j];
                            fzij +=  qi * T2.abD2[j][2] * muj[j];
                        }
                    }


                }
                /*
                   ===========================================================  
                   ===========================================================  
                */
                /* TOTAL FORCES */
                fx_dir[ia] += -fxij;
                fy_dir[ia] += -fyij;
                fz_dir[ia] += -fzij;
                fx_dir[ja] +=  fxij;
                fy_dir[ja] +=  fyij;
                fz_dir[ja] +=  fzij;
                /* TOTAL STRESS TENSOR */
                ttau[0][0] += - (rij[0] * fxij + rij[0] * fxij)*0.5;
                ttau[1][1] += - (rij[1] * fyij + rij[1] * fyij)*0.5;
                ttau[2][2] += - (rij[2] * fzij + rij[2] * fzij)*0.5;
                ttau[0][1] += - (rij[0] * fyij + rij[1] * fxij)*0.5;
                ttau[0][2] += - (rij[0] * fzij + rij[2] * fxij)*0.5;
                ttau[1][2] += - (rij[1] * fzij + rij[2] * fyij)*0.5;
            } /* ja */
        } /* ia */
    } /* omp */ 


    *u_dir = uu; 

    ttau[1][0]=ttau[0][1];
    ttau[2][0]=ttau[0][2];
    ttau[2][1]=ttau[1][2];
    for(int i=0; i<3;i++){
        for(int j=0; j<3;j++){
            tau_dir[i][j] =  ttau[i][j] * simuCell.inveOmegaPU;
        }
    }
    /*
    printf("tau_dir\n");
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            printf(ee,tau_dir[i][j]);
        }
        putchar('\n');
    }
    putchar('\n');
    */
//    for(int ia = 0; ia < nion; ia ++){
//        ef_dir[ia][0]=ef_x[ia];
//        ef_dir[ia][1]=ef_y[ia];
//        ef_dir[ia][2]=ef_z[ia];
//        for(int i=0; i<3;i++){
//            ef_dir[ia][i]=tef[ia][i];
//            for(int j=0; j<3;j++){
//                efg_dir[ia][i][j]=tefg[ia][i][j];
//            }
//        }
//    }

//    printf("Udmp %e\n",uudmp);
    /*************************************** 
            direct to cartesian                   
     ***************************************/
    dirkar ( nion , rx , ry , rz , simuCell.A ) ;

//    free(ef_x);
//    free(ef_y);
//    free(ef_z);
//    free(tefg);
}

/******************************************************************************/
/* Multipole expansion of the Ewald sum  in reciprocal space                  */
/******************************************************************************/
void multipole_ES_rec(double *q, double (*mu)[3], double (*theta)[3][3],
                      double *u_rec  , double (*ef_rec)[3], double (*efg_rec)[3][3],
                      double *fx_rec, double *fy_rec, double *fz_rec , double tau_rec[3][3],
                      bool lqchtask, bool ldiptask){

    double uu;
    double kx,ky,kz,Ak,kcoe;
    double rxi,ryi,rzi;
    double fxij,fyij,fzij;
    double qi;
    double k_dot_r,k_dot_mu;
    double ckria,skria;
    double rhonk_R,rhonk_I,str;
    double recarg;
    double recarg2;

    double ttau[3][3];
    for (int i=0;i<3;i++){
        for(int j=i;j<3;j++){
            ttau[i][j]=0.0;
        }
    }
    int ik,ia;

    uu=0;

    #pragma omp parallel shared(rx,ry,rz) \
                         private(ik,qi,kx,ky,kz,Ak,kcoe,rxi,ryi,rzi,k_dot_r,k_dot_mu,rhonk_R,rhonk_I,str,fxij,fyij,fzij,ckria,skria)
    { 
        #pragma omp for reduction (+:uu,ttau,fx_rec[:nion],fy_rec[:nion],fz_rec[:nion],ef_rec[:nion],efg_rec[:nion]) schedule(dynamic,8) 
        for (ik=kcoul.kptDec.iaStart;ik<kcoul.kptDec.iaEnd;ik++){
            //if (kcoul.kk[ik] == 0.0) continue;
            kx   = kcoul.kx[ik];
            ky   = kcoul.ky[ik];
            kz   = kcoul.kz[ik];
            Ak   = kcoul.Ak[ik];
            kcoe = kcoul.kcoe[ik];
    
            rhonk_R = 0.0;
            rhonk_I = 0.0;
    
            for (ia=0;ia<nion;ia++){
                qi  = q[ia];
                rxi = rx[ia];
                ryi = ry[ia];
                rzi = rz[ia];
                k_dot_r  = ( kx * rxi + ky * ryi + kz * rzi );
                ckria  = cos(k_dot_r);
                skria  = sin(k_dot_r);
                rhonk_R  += qi * ckria;
                rhonk_I  += qi * skria;
        
         //       if ( ! ldiptask ) continue;
                k_dot_mu = ( mu[ia][0] * kx + mu[ia][1] * ky + mu[ia][2] * kz );
                rhonk_R += -k_dot_mu * skria;
                rhonk_I +=  k_dot_mu * ckria;
        
                
            } /* sum to get charge density */

            str =  (rhonk_R*rhonk_R + rhonk_I*rhonk_I) * Ak;
            /* potential energy */                                       
            uu  += str ;                                  
        
            /* second sum */
            for (ia=0;ia<nion;ia++){
                qi=q[ia];
                rxi = rx[ia];
                ryi = ry[ia];
                rzi = rz[ia];
                k_dot_r  = ( kx * rxi + ky * ryi + kz * rzi );
                ckria  = cos(k_dot_r);
                skria  = sin(k_dot_r);
                recarg  = Ak * (rhonk_I*ckria - rhonk_R*skria);
            
                fxij = kx * recarg;
                fyij = ky * recarg;
                fzij = kz * recarg;
    
                ef_rec[ia][0] += -fxij ;
                ef_rec[ia][1] += -fyij ;
                ef_rec[ia][2] += -fzij ;
    
                recarg2 = Ak * ( rhonk_R * ckria + rhonk_I * skria );
                efg_rec[ia][0][0] += kx * kx * recarg2; 
                efg_rec[ia][1][1] += ky * ky * recarg2; 
                efg_rec[ia][2][2] += kz * kz * recarg2; 
                efg_rec[ia][0][1] += kx * ky * recarg2; 
                efg_rec[ia][0][2] += kx * kz * recarg2; 
                efg_rec[ia][1][2] += ky * kz * recarg2; 
    
                fx_rec[ia] += -qi * fxij;
                fy_rec[ia] += -qi * fyij;
                fz_rec[ia] += -qi * fzij;
        //        if (! ldiptask ) continue ;
                k_dot_mu = ( mu[ia][0] * kx + mu[ia][1] * ky + mu[ia][2] * kz ) * recarg2;
                fx_rec[ia] += kx * k_dot_mu;
                fy_rec[ia] += ky * k_dot_mu;
                fz_rec[ia] += kz * k_dot_mu;
    
    
            }
            ttau[0][0] +=   ( 1.0 - kcoe * kx * kx ) * str;
            ttau[0][1] +=         - kcoe * kx * ky   * str;
            ttau[0][2] +=         - kcoe * kx * kz   * str;
            ttau[1][1] +=   ( 1.0 - kcoe * ky * ky ) * str;
            ttau[1][2] +=         - kcoe * ky * kz   * str;
            ttau[2][2] +=   ( 1.0 - kcoe * kz * kz ) * str;
        }
    } /* omp parallel */

    //sample_("f_rec",fx_rec,fy_rec,fz_rec);

    /* half mesh and units */
    double ttpiV = 2.0 * TPI * simuCell.inveOmega;
    double tfpiV = 2.0 * ttpiV;
    uu = uu * ttpiV;
    *u_rec =   uu ;

    for (int ia=0;ia<nion;ia++){
        fx_rec[ia] *= tfpiV;
        fy_rec[ia] *= tfpiV;
        fz_rec[ia] *= tfpiV;
        for (int j=0;j<3;j++){
            ef_rec[ia][j] *= tfpiV;
            for (int k=0;k<3;k++){
                efg_rec[ia][j][k] *= tfpiV;
            }
        }
    }
    ttau[1][0]=ttau[0][1];
    ttau[2][0]=ttau[0][2];
    ttau[2][1]=ttau[1][2];
    for (int j=0;j<3;j++){
        for (int k=0;k<3;k++){
            tau_rec[j][k] = ttau[j][k] * ttpiV * simuCell.inveOmegaPU;
        }
    }
    /*
    printf("tau_rec\n");
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            printf(ee,tau_rec[i][j]);
        }
        putchar('\n');
    }
    putchar('\n');
    */

    //printf("at the end EFG %e %e\n",efg_rec[0][0][0],efg_rec[0][1][1]);
    /* remark on unit
    1/(4*pi*epislon_0) = 1 => epsilon_0 = 1/4pi
    */



}

