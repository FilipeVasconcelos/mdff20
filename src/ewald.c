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

#ifdef DEBUG
    #define DEBUG_EWALD
    #define DEBUG_EWALD_DIR
    #define DEBUG_EWALD_REC
    #define DEBUG_EWALD_DIR_DIPOLE
    #define DEBUG_EWALD_DIR_QUADRIPOLE
    #define DEBUG_EWALD_DIR_COUNT
    #define DEBUG_EWALD_REC_COUNT
    #define DEBUG_EWALD_DIR_PIM_COUNT
    #define DEBUG_EWALD_REC_PIM_COUNT
#endif
//#define DEBUG_
#ifdef DEBUG_
    #define DEBUG_EWALD
    #define DEBUG_EWALD_DIR_DIPOLE
    #define DEBUG_EWALD_DIR_QUADRIPOLE
    #define DEBUG_EWALD_DIR_COUNT
    #define DEBUG_EWALD_REC_COUNT
    #define DEBUG_EWALD_DIR_PIM_COUNT
    #define DEBUG_EWALD_REC_PIM_COUNT
#endif
//    #define DEBUG_EWALD_DIR
//    #define DEBUG_EWALD_REC

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
/* Ewald Summation */
void multipole_ES(double *q, double (*mu)[3], double (*theta)[3][3],double *u, double *pvir, double tau[3][3],
                             double (*ef)[3], double (*efg)[3][3], bool lq, bool ld, bool lt, bool lcouldamp,
                             bool do_forces, bool do_stress, bool do_ef, bool do_efg , bool do_dir, bool do_rec , 
                             bool inpim, bool update_sf ){

#ifdef DEBUG_EWALD
    printf("inside multipole_ES\n");
#endif
    double u_dir, u_rec;
    double (*ef_dir)[3], (*ef_rec)[3], (*ef_self)[3];
    double (*efg_dir)[3][3], (*efg_rec)[3][3], (*efg_self)[3][3];
    double *fx_dir, *fy_dir, *fz_dir;
    double *fx_rec, *fy_rec, *fz_rec;
    double tau_dir[3][3],tau_rec[3][3];


    if ( do_dir ) {
        statime(18);
        u_dir=0.0;
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
        multipole_ES_dir(q, mu, theta, &u_dir, ef_dir, efg_dir, fx_dir, fy_dir, fz_dir, tau_dir,
                         lq, ld, lt, lcouldamp,
                         do_forces, do_stress, do_ef, do_efg);

#ifdef DEBUG_EWALD_DIR
        for (int i=0;i<3;i++){
            printf("ef_dir  %f\n",ef_dir[0][i]);
            printf("ef_dir  %f\n",ef_dir[1][i]);
            printf("\n");
        }
#endif
        statime(19);
        if ( inpim ) {
            cc_multipole_ES_dir_pim+=1;
            mestime(&ewaldDirpimCPUtime,19,18);
        }
        else {
            cc_multipole_ES_dir+=1;
            mestime(&ewaldDirCPUtime,19,18);
        }
    } /* do DIR */

    /*************************************************************/
    /*                RECIPROCAL SPACE                           */
    /*************************************************************/
    if ( do_rec ) {
        statime(20);
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

        u_rec=0.0;
        multipole_ES_rec(q, mu, theta, &u_rec, ef_rec, efg_rec, fx_rec, fy_rec, fz_rec, tau_rec,
                         lq, ld, do_forces, do_stress, do_ef, do_ef, update_sf);

        statime(21);
        if ( inpim ) {
            cc_multipole_ES_rec_pim+=1;
            mestime(&ewaldRecpimCPUtime,21,20);
        }
        else {
            cc_multipole_ES_rec+=1;
            mestime(&ewaldRecCPUtime,21,20);
        }

#ifdef DEBUG_EWALD_REC
        for (int i=0;i<3;i++){
            printf("ef_rec  %f\n",ef_rec[0][i]);
            printf("ef_rec  %f\n",ef_rec[1][i]);
            printf("\n");
        }
#endif
    } /* do REC */

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
            ef[ia][i] = ef_self[ia][i];
#ifdef DEBUG_EWALD
            printf("EF : %d %d self : %15.8e ",ia,i,ef_self[ia][i]);
#endif
            if ( do_dir ) {
                ef[ia][i] += ef_dir[ia][i];
#ifdef DEBUG_EWALD
                printf("dir: %15.8e ",ef_dir[ia][i]);
#endif
            }
            if ( do_rec ) {
                ef[ia][i] += ef_rec[ia][i];
#ifdef DEBUG_EWALD
                printf("rec: %15.8e ",ef_rec[ia][i]);
#endif
            }
#ifdef DEBUG_EWALD
                printf("tot: %15.8e \n",ef[ia][i]);
#endif
        /* ------------------------- */
        /* electric field gradient   */
        /* ------------------------- */
            for (int j=0;j<3;j++){
                efg[ia][i][j]=efg_self[ia][i][j];
                if ( do_dir ) {
                    efg[ia][i][j] += efg_dir[ia][i][j];
                }
                if ( do_rec ) {
                    efg[ia][i][j] += efg_rec[ia][i][j];
                }
            }
        }
        /* ---------------------------- */
        /*           forces             */
        /* ---------------------------- */

        if ( do_dir ) {
            fx[ia] += fx_dir[ia] * coul_unit;
            fy[ia] += fy_dir[ia] * coul_unit;
            fz[ia] += fz_dir[ia] * coul_unit;
        }
        if ( do_rec ) {
            fx[ia] += fx_rec[ia] * coul_unit;
            fy[ia] += fy_rec[ia] * coul_unit;
            fz[ia] += fz_rec[ia] * coul_unit;
        }
    }

    /* ---------------------------- */
    /*         stress tensor        */
    /* ---------------------------- */
    if ( do_stress ) {
        for (int i=0;i<3;i++){
            for (int j=0;j<3;j++){
                tau[i][j]=0.0;
                if ( do_dir ) {
                    tau[i][j] += tau_dir[i][j] * coul_unit;
                }
                if ( do_rec ) {
                    tau[i][j] += tau_rec[i][j] * coul_unit;
                }
            }
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
    printf("u_self %e\n",u_self*coul_unit);
    printf("u_coul %e\n",*u);
    if ( do_dir ) {
        printf("u_dir  %e\n",u_dir*coul_unit);
        printf("fx_dir  %e %e %e \n",fx_dir[0],fy_dir[0],fz_dir[0]);
        printf("ef_dir  %e %e %e \n",ef_dir[0][0],ef_dir[0][1],ef_dir[0][2]);
    }
    if ( do_rec ) {
        printf("u_rec  %e\n",u_rec*coul_unit);
        printf("fx_rec  %e %e %e \n",fx_rec[0],fx_rec[0],fx_rec[0]);
        printf("ef_rec  %e %e %e\n",ef_rec[0][0],ef_rec[0][1],ef_rec[0][2]);
        putchar('\n');
    }
#endif

    free(ef_self);
    free(efg_self);
    if ( do_dir ) {
        free(ef_dir);
        free(fx_dir);
        free(fy_dir);
        free(fz_dir);
        free(efg_dir);
    }
    if ( do_rec ) {
        free(ef_rec);
        free(efg_rec);
        free(fx_rec);
        free(fy_rec);
        free(fz_rec);
    }
#ifdef DEBUG_EWALD
    printf("outside multipole_ES\n");
#endif
}


/******************************************************************************/
/* Multipole expansion of the Ewald sum  in direct space                      */
void multipole_ES_dir(double *q, double (*mu)[3], double (*theta)[3][3],
                      double *u_dir, double (*ef_dir)[3], double (*efg_dir)[3][3],
                      double *fx_dir, double *fy_dir, double *fz_dir , double tau_dir[3][3],
                      bool lqchtask, bool ldiptask, bool lquatask, bool lcouldamp,
                      bool do_forces, bool do_stress, bool do_ef, bool do_efg){

#ifdef DEBUG_EWALD_DIR
    printf("inside multipole_ES_dir\n");
#endif
    double rxi,ryi,rzi;
    double rij[3];
    double fij[3];
    TENSOR_RK0 T0;
    TENSOR_RK1 T1;
    TENSOR_RK2 T2;
    TENSOR_RK3 T3;
    TENSOR_RK4 T4;
    TENSOR_RK5 T5;
    double qi;
    double mui[3];
    double thetai[3][3];
    double qj;
    double muj[3];
    bool ldamp;
    bool qch_i,qch_j,qch_iOUj;
    bool dip_i,dip_j,dip_iETj,dip_iOUj;
    bool qua_i,qua_j,qua_iOUj,qua_iETj;
    double thetaj[3][3];

    double expon;
    double qij;
    double d , d2 , d3  , d5 , d7, d9;
    double dm1 , dm3 , dm5  , dm7 , dm9, dm11;
    double F0, F1, F2, F3, F4, F5;
    double F1_dm3 ,F2_dm5 , F3_dm7 , F4_dm9, F5_dm11;

    // damping related
    double fdamp1,fdamp2,fdampdiff1,fdampdiff2;
    double F1d1, F2d1, F1d2, F2d2;
    double F2d1_dm5,F2d2_dm5,F1d1_dm3,F1d2_dm3;

    double alpha2, alpha3, alpha5, alpha7, alpha9;
    double uu = 0;
    int ia,j1,ja,jb,je,ita,jta;
    double ttau[3][3];
    for (int i=0;i<3;i++){
        mui[i]=0.0;
        muj[i]=0.0;
        for(int j=i;j<3;j++){
            ttau[i][j]=0.0;
            thetai[i][j]=0.0;
            thetaj[i][j]=0.0;
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


    #pragma omp parallel default(none) \
    shared(lrcutsq,q,mu,theta,rx,ry,rz,typia,piroot,pol_damp_b,pol_damp_c,pol_damp_k,lpoldamping,nion,\
          efg_dir,ef_dir,lqchdiptask,lqchtask,ldiptask,lquatask,alphaES,alpha2,alpha3,alpha5,alpha7,alpha9,\
          verlet_coul,fx_dir,fy_dir,fz_dir,ttau,uu,atomDec,lverletL,lcouldamp,atypia,do_forces,do_stress,do_ef,do_efg) \
    private(ita,jta,ia,j1,jb,je,ja,qi,qj,qij,mui,muj,thetai,thetaj,rij,rxi,ryi,rzi,d2,qch_i,dip_i,qua_i,qch_j,dip_j,qua_j,\
            dip_iETj,dip_iOUj,qch_iOUj,qua_iOUj,qua_iETj,d,d3,d5,d7,d9,dm1,dm3,dm5,dm7,dm9,dm11,T0,T1,T2,T3,T4,T5,\
            fdamp1,fdamp2,fdampdiff1,fdampdiff2,expon,F0,F1,F2,F3,F4,F5,F1d1,F2d1,F1d2,F2d2,\
            F1_dm3,F2_dm5,F3_dm7,F4_dm9,F5_dm11,ldamp,F2d1_dm5,F2d2_dm5,F1d1_dm3,F1d2_dm3,fij)
    {
        #pragma omp for reduction (+:uu,ttau,fx_dir[:nion],fy_dir[:nion],fz_dir[:nion],\
                                   ef_dir[:nion],efg_dir[:nion]) schedule(dynamic,8)
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
                for(int j=0;j<3;j++){
                    thetai[i][j] = theta[ia][i][j];
                    if (thetai[i][j] != 0.0) qua_i=true;
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
                dip_j=false;
                jta=typia[ja];
                for(int i=0;i<3;i++){
                    fij[i] = 0.0;
                    muj[i] = mu[ja][i];
                    if (muj[i] != 0.0) dip_j=true;
                    for(int j=0;j<3;j++){
                        thetaj[i][j] = theta[ja][i][j];
                        if (thetaj[i][j] != 0.0) qua_j=true;
                    }
                }
                qij    = qi * qj ;
                rij[0] = rxi - rx[ja];
                rij[1] = ryi - ry[ja];
                rij[2] = rzi - rz[ja];
                pbc(&rij[0],&rij[1],&rij[2]);
                d2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
                if ( d2 > lrcutsq ) continue ;
                qch_iOUj = qch_i || qch_j;
                dip_iOUj = dip_i || dip_j;
                dip_iETj = dip_i && dip_j;
                qua_iOUj = qua_i || qua_j;
                qua_iETj = qua_i && qua_j;
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

		// damping function
                ldamp = ( (lpoldamping[ita][ita][jta] || lpoldamping[jta][ita][jta] ) && ( lcouldamp ) );
                if ( ldamp ) {
#ifdef DEBUG_DIPOLE_DAMPING
        printf("damping: %d %d %d %d %d %e %e %d\n",ia,ja,ita,jta,lcouldamp,pol_damp_b[ita][ita][jta],pol_damp_c[ita][ita][jta],pol_damp_k[ita][ita][jta]);
#endif
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
            		        for(int m=0; m<3;m++){
                                    T5.abcde[i][j][k][l][m]=0.0;
                                }
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
                    T1.a[i] = - rij[i] * dm3;
                    if ( ldamp ) {
                        T1.aD1[i] = T1.a[i] * F1d1;
                        T1.aD2[i] = T1.a[i] * F1d2;
                    }
                    T1.a[i] *= F1;
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
                if (lquatask) {
                    /*******************************************/
                    /* multipole interaction tensor rank = 5   */
                    /* nb of components = 243 => reduced = ?   */
                    /*******************************************/
                    F5_dm11 = dm11 * F5 * 945.0;
    	            for(int i=0; i<3;i++){
                        for(int j=0; j<3;j++){
                            for(int k=0; k<3;k++){
                                for(int l=0; l<3;l++){
                                    for(int m=0; m<3;m++){
                                        T5.abcde[i][j][k][l][m] = rij[i] * rij[j] * rij[k] * rij[l] * rij[m] * F5_dm11;
                                        if ( l == m ) T5.abcde[i][j][k][l][m] += - rij[i]*rij[j]*rij[k] * F4_dm9;
                                        if ( k == m ) T5.abcde[i][j][k][l][m] += - rij[i]*rij[j]*rij[l] * F4_dm9;
                                        if ( k == l ) T5.abcde[i][j][k][l][m] += - rij[i]*rij[j]*rij[m] * F4_dm9;
                                        if ( j == m ) T5.abcde[i][j][k][l][m] += - rij[i]*rij[k]*rij[l] * F4_dm9;
                                        if ( j == l ) T5.abcde[i][j][k][l][m] += - rij[i]*rij[k]*rij[m] * F4_dm9;
                                        if ( j == k ) T5.abcde[i][j][k][l][m] += - rij[i]*rij[l]*rij[m] * F4_dm9;
                                        if ( i == m ) T5.abcde[i][j][k][l][m] += - rij[j]*rij[k]*rij[l] * F4_dm9;
                                        if ( i == l ) T5.abcde[i][j][k][l][m] += - rij[j]*rij[k]*rij[m] * F4_dm9;
                                        if ( i == k ) T5.abcde[i][j][k][l][m] += - rij[j]*rij[l]*rij[m] * F4_dm9;
                                        if ( i == j ) T5.abcde[i][j][k][l][m] += - rij[k]*rij[l]*rij[m] * F4_dm9;
                                        if ( ( i == j ) && ( k == l ) ) T5.abcde[i][j][k][l][m] += rij[m] * F3_dm7;
                                        if ( ( i == j ) && ( l == m ) ) T5.abcde[i][j][k][l][m] += rij[k] * F3_dm7;
                                        if ( ( i == j ) && ( k == m ) ) T5.abcde[i][j][k][l][m] += rij[l] * F3_dm7;
                                        if ( ( i == k ) && ( j == l ) ) T5.abcde[i][j][k][l][m] += rij[m] * F3_dm7;
                                        if ( ( i == k ) && ( l == m ) ) T5.abcde[i][j][k][l][m] += rij[j] * F3_dm7;
                                        if ( ( i == k ) && ( j == m ) ) T5.abcde[i][j][k][l][m] += rij[l] * F3_dm7;
                                        if ( ( i == l ) && ( k == m ) ) T5.abcde[i][j][k][l][m] += rij[j] * F3_dm7;
                                        if ( ( i == l ) && ( j == m ) ) T5.abcde[i][j][k][l][m] += rij[k] * F3_dm7;
                                        if ( ( j == k ) && ( l == m ) ) T5.abcde[i][j][k][l][m] += rij[i] * F3_dm7;
                                        if ( ( j == k ) && ( i == m ) ) T5.abcde[i][j][k][l][m] += rij[l] * F3_dm7;
                                        if ( ( j == k ) && ( i == l ) ) T5.abcde[i][j][k][l][m] += rij[m] * F3_dm7;
                                        if ( ( j == l ) && ( k == m ) ) T5.abcde[i][j][k][l][m] += rij[i] * F3_dm7;
                                        if ( ( j == l ) && ( i == m ) ) T5.abcde[i][j][k][l][m] += rij[k] * F3_dm7;
                                        if ( ( k == l ) && ( j == m ) ) T5.abcde[i][j][k][l][m] += rij[i] * F3_dm7;
                                        if ( ( k == l ) && ( i == m ) ) T5.abcde[i][j][k][l][m] += rij[j] * F3_dm7;
                                    }
                                }
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
                if ( (lqchtask) && (qch_iOUj) ) {
#ifdef DEBUG_EWALD_DIR_
                    printf("inside charge-charge interaction %d %d %d %d\n",ia,ja,atomDec.iaStart,atomDec.iaEnd);
#endif
                    /*  potential energy */
                    uu += qij * T0.sca;
                    for(int i=0;i<3;i++){
                        // electric field
                        if ( do_ef ) {
                	    ef_dir[ia][i] += -qj * T1.a[i];
                            ef_dir[ja][i] +=  qi * T1.a[i];
                            if ( ldamp ) {
                                ef_dir[ia][i] +=  qj * T1.aD1[i];
                                ef_dir[ja][i] += -qi * T1.aD2[i];
                            }
                        }
                        /* forces */
                        if ( do_forces ){
                            fij[i] = qij * T1.a[i];
                        }
                        // electric field gradient
                        if ( do_efg ) {
                            for(int j=0;j<3;j++){
                                efg_dir[ia][i][j] += -qj * T2.ab[i][j];
                                efg_dir[ja][i][j] += -qi * T2.ab[i][j];
                            }
                        }
                    }
                }
                /*
                   ===========================================================
                                    dipole-dipole interaction
                   ===========================================================
                */
                if ( ( ldiptask ) && (dip_iETj) ) {
#ifdef DEBUG_EWALD_DIR_
                    printf("inside dipole-dipole interaction  %d %d\n",ia,ja);
#endif
                    for(int i=0; i<3;i++){
                        for(int j=0; j<3;j++){
                            /*  potential energy */
                            uu += - mui[i] * T2.ab[i][j] * muj[j];
                            /* electric field */
                            if ( do_ef ) {
                	            ef_dir[ia][i] += muj[j] * T2.ab[i][j] ;
                                ef_dir[ja][i] += mui[j] * T2.ab[i][j] ;
                            }
                            if ( (do_efg) || (do_forces) ) {
                                for(int k=0; k<3;k++){
                                    /* electric field gradient */
                                    efg_dir[ia][i][j] +=  muj[k] * T3.abc[i][j][k];
                                    efg_dir[ja][i][j] += -mui[k] * T3.abc[i][j][k];
                                    /* forces */
                                    fij[k] += -mui[i] * T3.abc[i][k][j] * muj[j];
                                }
                            }
                        }
                    }
                }/* ldiptask */
                /*
                   ===========================================================
                                    charge-dipole interaction
                   ===========================================================
                */
                if ( ( lqchdiptask ) && (dip_iOUj) ) {
#ifdef DEBUG_EWALD_DIR_
                    printf("inside charge-dipole interaction %d %d\n",ia,ja);
#endif
                    for(int i=0; i<3;i++){
                        /*  potential energy */
                        uu += - qi * T1.a[i] * muj[i];
                        uu +=   qj * T1.a[i] * mui[i];
                        /*  forces
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
                        if ( do_forces ) {
                            for(int j=0; j<3;j++){
                                fij[j] += -qi * T2.ab[j][i] * muj[i];
                                fij[j] +=  qj * T2.ab[j][i] * mui[i];
                            }
                        }
                        /* damping */
                        if (ldamp) {
                            /*  potential energy */
                            uu +=  qi *  T1.aD2[i] * muj[i];
                            uu += -qj *  T1.aD1[i] * mui[i];
                            if (do_forces) {
                                for(int j=0; j<3;j++){
                                    /*  forces */
                                    fij[j] += -qj * T2.abD1[i][j] * mui[i];
                                    fij[j] +=  qi * T2.abD2[i][j] * muj[i];
                                }
                            }
                        }
                    }

                } /* charge-dipole */

                /* TOTAL FORCES */
                if (do_forces) {
                    fx_dir[ia] += -fij[0];
                    fy_dir[ia] += -fij[1];
                    fz_dir[ia] += -fij[2];
                    fx_dir[ja] +=  fij[0];
                    fy_dir[ja] +=  fij[1];
                    fz_dir[ja] +=  fij[2];
                }
                /* TOTAL STRESS TENSOR */
                if (do_stress) {
                    for(int i=0;i<3;i++){
                        for(int j=i;j<3;j++){
                            ttau[i][j] += - (rij[i]*fij[j] + rij[j]*fij[i]);
                        }
                    }
                }
            } /* ja */
        } /* ia */
    } /* omp */


    *u_dir = uu;

    if ( do_stress ) {
        ttau[1][0]=ttau[0][1];
        ttau[2][0]=ttau[0][2];
        ttau[2][1]=ttau[1][2];
        for(int i=0; i<3;i++){
            for(int j=0; j<3;j++){
                tau_dir[i][j] =  0.5*ttau[i][j] * simuCell.inveOmegaPU;
            }
        }
    }

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
                      double *fx_rec, double *fy_rec, double *fz_rec , double tau_rec[3][3],
                      bool lqchtask, bool ldiptask,bool do_forces, bool do_stress, bool do_ef, bool do_efg,
                      bool update_sf){

#ifdef DEBUG_EWALD_REC
    printf("inside multipole_ES_rec\n");
#endif
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
    if ( (update_sf) && (lpim) ) struct_fact();

    #pragma omp parallel default(none) \
                         shared(rx,ry,rz,q,mu,theta,u_rec,ef_rec,efg_rec,fx_rec,fy_rec,fz_rec,tau_rec,\
                                lqchtask,ldiptask,do_ef,do_efg,do_forces,do_stress,kcoul,uu,ttau,nion,lpim) \
                         private(ik,ia,qi,kx,ky,kz,Ak,kcoe,rxi,ryi,rzi,k_dot_r,k_dot_mu,rhonk_R,rhonk_I,str,fxij,fyij,fzij,ckria,skria,recarg,recarg2)
    {
        #pragma omp for reduction (+:uu,ttau,fx_rec[:nion],fy_rec[:nion],fz_rec[:nion],ef_rec[:nion],efg_rec[:nion]) \
                        schedule(dynamic,8)
        for (ik=kcoul.kptDec.iaStart;ik<kcoul.kptDec.iaEnd;ik++){

            if (kcoul.kk[ik] == 0.0) continue;
            kx   = kcoul.kx[ik];
            ky   = kcoul.ky[ik];
            kz   = kcoul.kz[ik];
            Ak   = kcoul.Ak[ik];
            kcoe = kcoul.kcoe[ik];
            kcoul.rhon_R[ik]=0.0;
            kcoul.rhon_I[ik]=0.0;

            if ( lpim ) { 
                /* charge density at ik */ 
                if ( ( lqchtask ) && ( ldiptask ) ) {
                    charge_density_qmu(ik,q,mu);
                } else {
                    if ( lqchtask ) charge_density_q(ik,q);
                    if ( ldiptask ) charge_density_mu(ik,mu);
                }
                rhonk_R = kcoul.rhon_R[ik];
                rhonk_I = kcoul.rhon_I[ik];
            } 
            else {
                struct_fact_rhon(ik,q,mu,lqchtask,ldiptask,&rhonk_R,&rhonk_I);
                //rhonk_R = kcoul.rhon_R[ik];
                //rhonk_I = kcoul.rhon_I[ik];
            }
#ifdef DEBUG_EWALD_REC
            printf("k %d %d %d"ee3"\n",ik,kcoul.kptDec.iaStart,kcoul.kptDec.iaEnd,str,rhonk_R,rhonk_I);
#endif
            str =  (rhonk_R*rhonk_R + rhonk_I*rhonk_I) * Ak;
            /* potential energy */
            uu  += str ;

            /* second sum */
            for (ia=0;ia<nion;ia++){
                qi=q[ia];
                ckria = kcoul.ckria[ik][ia];
                skria = kcoul.skria[ik][ia];
                /*
                rxi = rx[ia];
                ryi = ry[ia];
                rzi = rz[ia];
                k_dot_r  = ( kx * rxi + ky * ryi + kz * rzi );
                ckria=cos(k_dot_r);
                skria=sin(k_dot_r);
                */

                if ( ( do_ef ) || ( do_forces ) ) {
                    recarg  = Ak * (rhonk_I*ckria - rhonk_R*skria);
                    fxij = kx * recarg;
                    fyij = ky * recarg;
                    fzij = kz * recarg;
                }

                if ( do_ef ) {
                    ef_rec[ia][0] += -fxij ;
                    ef_rec[ia][1] += -fyij ;
                    ef_rec[ia][2] += -fzij ;
                }

                if ( do_efg ) {
                    recarg2 = Ak * ( rhonk_R * ckria + rhonk_I * skria );
                    efg_rec[ia][0][0] += kx * kx * recarg2;
                    efg_rec[ia][1][1] += ky * ky * recarg2;
                    efg_rec[ia][2][2] += kz * kz * recarg2;
                    efg_rec[ia][0][1] += kx * ky * recarg2;
                    efg_rec[ia][0][2] += kx * kz * recarg2;
                    efg_rec[ia][1][2] += ky * kz * recarg2;
                }

                if ( do_forces ) {
                    fx_rec[ia] += -qi * fxij;
                    fy_rec[ia] += -qi * fyij;
                    fz_rec[ia] += -qi * fzij;
                    if (! ldiptask ) continue ;
                    k_dot_mu = ( mu[ia][0] * kx + mu[ia][1] * ky + mu[ia][2] * kz ) * recarg2;
                    fx_rec[ia] += kx * k_dot_mu;
                    fy_rec[ia] += ky * k_dot_mu;
                    fz_rec[ia] += kz * k_dot_mu;
                }


            }
            if ( do_stress ) {
                ttau[0][0] +=   ( 1.0 - kcoe * kx * kx ) * str;
                ttau[0][1] +=         - kcoe * kx * ky   * str;
                ttau[0][2] +=         - kcoe * kx * kz   * str;
                ttau[1][1] +=   ( 1.0 - kcoe * ky * ky ) * str;
                ttau[1][2] +=         - kcoe * ky * kz   * str;
                ttau[2][2] +=   ( 1.0 - kcoe * kz * kz ) * str;
            }
        }
    } /* omp parallel */

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

}

