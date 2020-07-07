#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "io.h"
#include "tools.h"
#include "field.h"
#include "pim.h"
#include "config.h"
#include "md.h"
#include "ewald.h"
#include "coulombic.h"

#ifdef DEBUG
    #define DEBUG_ASPC
    #define DEBUG_INNER_KO
#endif
//#define DEBUG_
#ifdef DEBUG_
    #define DEBUG_ASPC
    #define DEBUG_INNER_KO
#endif
/******************************************************************************/
int read_pim (char * controlfn)
{
    char buffer[MAX_LEN+1];
    FILE * fp;

    fp = fopen (controlfn, "r");

    if (NULL == fp ) {
        pError("opening control file (reading pim) ");
        return (-1);
    }

    while (EOF != fscanf(fp, "%s\n", buffer)) {
        if (strcmp(buffer,"algo_pim") == 0 ) {
            fscanf(fp,"%s",buffer);
            algo_pim=check_string("algo_pim",buffer,allwd_algo_pim, ALLWD_ALGO_PIM_STR);
        }
        if (strcmp(buffer,"algo_extrapolate_dipole") == 0 ) {
            fscanf(fp,"%s",buffer);
            algo_extrapolate_dipole=check_string("algo_extrapolate_dipole",buffer,allwd_algo_extrapolate_dipole,
                                                           ALLWD_ALGO_EXTRAPOLATE_DIPOLE_STR);
        }
        if (strcmp(buffer,"conv_tol_ind") == 0 ) {
            fscanf(fp,"%lf",&conv_tol_ind);
        }
        if (strcmp(buffer,"extrapolate_order") == 0 ) {
            fscanf(fp,"%d",&extrapolate_order);
        }
        if (strcmp(buffer,"min_scf_pim_iter") == 0 ) {
            fscanf(fp,"%d",&min_scf_pim_iter);
        }
        if (strcmp(buffer,"max_scf_pim_iter") == 0 ) {
            fscanf(fp,"%d",&max_scf_pim_iter);
        }
        if (strcmp(buffer,"omegakO") == 0 ) {
            fscanf(fp,"%lf",&omegakO);
        }
   }
   fclose(fp);
   return 0;
}

/******************************************************************************/
void default_pim(){
    /* gen allowed input strings for algoPIM */
    strcpy(allwd_algo_pim[0],"scf");
    strcpy(allwd_algo_pim[1],"scfkO");
    /* gen allowed input strings for algo_ext_dipole */
    strcpy(allwd_algo_extrapolate_dipole[0],"poly");
    strcpy(allwd_algo_extrapolate_dipole[1],"aspc");
    conv_tol_ind  = 1e-6;
    min_scf_pim_iter = 3;
    max_scf_pim_iter = 100;
    extrapolate_order = 0 ;
    algo_extrapolate_dipole = 0;
    algo_pim = 0;
}

/******************************************************************************/
void check_pim(){
}

/******************************************************************************/
void init_pim (char * controlfn) {

    default_pim();
    read_pim(controlfn);
    info_pim();
}

/******************************************************************************/
void info_pim(){
    if (ionode) {
        SEPARATOR;
        printf("pim info\n");
        LSEPARATOR;
        putchar('\n');
        printf("Polarizable Ion Model\n");
        switch (algo_pim) {
            case 0 :
                printf("SCF algorithm\n");
                putchar('\n');
                printf("minimum nb iterations : %d\n",min_scf_pim_iter);
                printf("maximum nb iterations : %d\n",max_scf_pim_iter);
                break;
            case 1 :
                printf("SCF algorithm by kO\n");
                printf("minimum nb iterations : %d\n",min_scf_pim_iter);
                printf("maximum nb iterations : %d\n",max_scf_pim_iter);
                printf("omegakO               : %e\n",omegakO);

                break;
        }
        switch (algo_extrapolate_dipole) {
            case 0 :
                printf("polynomial extrapolation\n");
                break;
            case 1 :
                printf("ASPC (Always Stable Predictor Corrector) extrapolation\n");
                putchar('\n');
                printf("extrapolate order : %d\n",extrapolate_order);
                break;
        }
        printf("convergence tolerance : %e\n",conv_tol_ind);
    }
}

/*****************************************************************************
     ---------------------------------------------------------------
          u_pol = 1/ ( 2 alpha  ) * | mu_ind | ^ 2
     ---------------------------------------------------------------
     note on units :
        [ poldipia  ] = A ^ 3
        [ mu_ind ] = e A
        [ u_pol  ] = e^2 / A ! electrostatic internal energy
*****************************************************************************/
double get_upol(double (*mu)[3]){
    double upol=0.0;
    for(int ia=0;ia<nion;ia++){
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                upol+= mu[ia][i]*invepolia[ia][i][j]*mu[ia][j];
#ifdef NON_DIAGONAL_POLIA
                if ( fabs(polia[ia][i][j]) > 1e-12 ) {
                   printf("test invepolia %e\n",1.0/polia[ia][i][j]);
                    upol+= mu[ia][i]*(1.0/polia[ia][i][j])*mu[ia][j];
                }
#endif
            }
        }
    }
    return upol*0.5;
}

/******************************************************************************/
/* simple implementation of the self consistent algorithm :                   */
/*
 * Main Goal is to get induced electric field self consistently               */
/* Step 1 : Get the electrostatic field (ef_stat)
 *          from static charges and dipoles                                   */
/* Step 2 : Get the induced electric dipole (mu_ind) (self consistent)
 *          from electric field ef (ef=ef_stat) for the first iteration       */
/* Step 3 : Get the induced electric field (ef_ind) from induced
 *          electric dipole only
 * Step 4 : Calculate ef = ef_stat + ef_ind                                   */
/* Step 5 : if not converged go to Step 2                                     */
/*                                                                            */
/******************************************************************************/
void momentpolaSCF(double (*mu_ind)[3],double *upol){

#ifdef DEBUG_SCF
    printf("inside momentpolaSCF \n");
#endif
    int iscf;
    double rmsd;
    double (*ef)[3];
    double (*ef_stat)[3];
    double (*ef_ind)[3];
    double u_coul_stat,u_coul_ind;
    double u_coul_pol;
    double uupol;
    double pdummy;
    double tdummy[3][3];
    double (*efgdummy)[3][3];
    double (*zeroch);
    zeroch=malloc(nion*sizeof(*zeroch));
    ef=malloc(nion*sizeof(*ef));
    efgdummy=malloc(nion*sizeof(*efgdummy));
    ef_stat=malloc(nion*sizeof(*ef_stat));
    ef_ind=malloc(nion*sizeof(*ef_ind));
    for(int ia=0;ia<nion;ia++){
        zeroch[ia]=0.0;
        for(int i=0;i<3;i++){
            ef[ia][i]=0.0;
            ef_stat[ia][i]=0.0;
            ef_ind[ia][i]=0.0;
            for(int j=0;j<3;j++){
                efgdummy[ia][i][j]=0.0;
                tdummy[i][j]=0.0;
            }
        }
    }

    /* static electric field from
     * charges+dipoles+quadrupoles
     * and damping
    */
    multipole_ES(qia,dipia,quadia,&u_coul_stat,&pdummy,tdummy,ef_stat,efgdummy,lqch,ldip,lqua,ldmp,
                /* do   forces stress ef    efg    dir   rec   inpim*/
                        false, false, true, false, true, true, true );

    /* first electric field is static */
    for (int ia=0;ia<nion;ia++){
        for (int i=0;i<3;i++){
            ef[ia][i]=ef_stat[ia][i];
        }
    }
    /* SCF LOOP */
    iscf=0;
    rmsd=DBL_MAX;
    if ( iopnode(istep,npas,nprint) ) {
        printf("  -----------------------------------------------------------\n");
        printf("                 PIM :     self consistent                   \n");
        printf("  -----------------------------------------------------------\n");
        printf("        iter          u_pol          u_coul            rmsd  \n");
        printf("  -----------------------------------------------------------\n");
    }

    while ( ( (iscf<max_scf_pim_iter) && (rmsd>conv_tol_ind))  || (iscf<min_scf_pim_iter) ) {

        for (int ia=0;ia<nion;ia++){
            mu_ind[ia][0]=0.0;
            mu_ind[ia][1]=0.0;
            mu_ind[ia][2]=0.0;
            ef_ind[ia][0]=0.0;
            ef_ind[ia][1]=0.0;
            ef_ind[ia][2]=0.0;
        }
        /* predictor */
        if ( ! iscf ) {
            extrapolate_dipole_aspc(mu_ind,ef,1);
        }
        else {
            induced_moment(mu_ind,ef);
        }
        uupol=get_upol(mu_ind);

        /* induced electric field from only induced dipoles */
        /* quadia -> quad_ind if induce electric field gradient is wanted  */
        /*                                                                   qch   dip  quad  dmp */
        multipole_ES(zeroch,mu_ind,quadia,&u_coul_ind,&pdummy,tdummy,ef_ind,efgdummy,false,true,false,false,
        /* do   forces stress ef    efg    dir   rec  inpim */
                false, false, true, false, true, true, true );

        /* ef = ef_stat + ef_ind */
        for (int ia=0;ia<nion;ia++){
            for (int i=0;i<3;i++){
                ef[ia][i] = ef_stat[ia][i] + ef_ind[ia][i];
            }
        }
        u_coul_pol = u_coul_stat - u_coul_ind;

        /* corrector */
        if (!iscf) extrapolate_dipole_aspc(mu_ind,ef,2);


        rmsd=get_rmsd_mu(mu_ind, ef);
        if ( iopnode(istep,npas,nprint ) )  {
            printf("  scf: %5d  "ee3"\n",iscf+1,uupol*coul_unit,u_coul_pol*coul_unit,rmsd);
        }
        iscf+=1;
    }

    //store induced dipole at t
    for (int ia=0;ia<nion;ia++){
        for (int i=0;i<3;i++){
            dipia_ind[ia][0][i]=mu_ind[ia][i];
        }
    }
    if ( iopnode(istep,npas,nprint ) )  {
        putchar('\n');
        printf("  scf converged in %d iterations, rmsd %e\n",iscf,rmsd);
        putchar('\n');
#ifdef DEBUG_ASPC
        printf("  ASPC scf dipia_ind ");
        for (int k=0;k<extrapolate_order;k++){
            printf("%e ",dipia_ind[0][k][0]);
        }
        printf("-> mu_ind %e\n",mu_ind[0][0]);
#endif
    }

    *upol=uupol*coul_unit;
    free(ef);
    free(efgdummy);
    free(ef_stat);
    free(ef_ind);
    free(zeroch);

}

/******************************************************************************/
/*  order of extrapolation is k+1 of Kolafa original derivation.              */
/*  zero order is simply the previous step                                    */
/*  if key == 1 => predictor                                                  */
/*  if key == 2 => corrector                                                  */
/******************************************************************************/
void extrapolate_dipole_aspc(double (*mu_ind)[3] , double (*ef)[3], int key ){

    double B_ASPC[MAX_EXTRAPOLATE_ORDER+1], W_ASPC;
    int ext_ord;

    double (*mu_save)[3];
    mu_save=malloc(nion*sizeof(*mu_save));

    // switch to lower order of extrapolation
    // if time step is not large enough to store previous dipoles steps
    if ( istep >= extrapolate_order+1 ) {
        ext_ord = extrapolate_order;
    }
    else {
        if ( ! istep ) { //if istep istep==1 or 0
            ext_ord = 0;
        }
        else {
          ext_ord = istep - 1;
        }
    }
    switch (ext_ord) {
        default:
            io_node pError("value of aspc extrapolation order not available should be 0<= extrapolate_order <= 4\n");
            exit(-1);
            break;
        case 0:
            B_ASPC[0] = 1.0;
            W_ASPC    = 0.0;
            break;
        case 1:
            B_ASPC[0] =  2.0;
            B_ASPC[1] = -1.0;
            W_ASPC    =  2.0/3.0;
            break;
        case 2:
            B_ASPC[0] =  2.5;
            B_ASPC[1] = -2.0;
            B_ASPC[2] =  0.5;
            W_ASPC    =  0.6;
            break;
        case 3:
            B_ASPC[0] =  2.8;
            B_ASPC[1] = -2.8;
            B_ASPC[2] =  1.2;
            B_ASPC[3] = -0.2;
            W_ASPC    =  4.0/7.0;
            break;
        case 4:
            B_ASPC[0] =  3.0;
            B_ASPC[1] = -24.0/7.0;
            B_ASPC[2] =  27.0/14.0;
            B_ASPC[3] = -4.0/7.0;
            B_ASPC[4] =  1.0/14.0;
            W_ASPC    =  5.0/9.0;
            break;
        case 5 :
            B_ASPC[0] =  22.0/7.0;
            B_ASPC[1] = -55.0/14.0;
            B_ASPC[2] =  55.0/21.0;
            B_ASPC[3] = -22.0/21.0;
            B_ASPC[4] =  5.0/21.0;
            B_ASPC[5] = -1.0/42.0;
            W_ASPC    =  6.0/11.0;
            break;
    }

    /* predictor */
    if (key == 1) {

        for (int ia=0;ia<nion;ia++){
            for (int i=0;i<3;i++){
                mu_ind[ia][i]=0.0;
            }
            for (int k=0;k<ext_ord+1;k++){
                for (int i=0;i<3;i++){
                    mu_ind[ia][i]+=B_ASPC[k]*dipia_ind[ia][k][i];
                }
            }
            if (ext_ord != extrapolate_order) ext_ord+=1;
            for (int k=ext_ord;k>0;k--){
                for (int i=0;i<3;i++){
                    dipia_ind[ia][k][i]=dipia_ind[ia][k-1][i];
                }
            }
        }
    }

    /* corrector */
    if (key == 2) {
        /* -------------------------*/
        /* save induced moment      */
        for (int ia=0;ia<nion;ia++){
            for (int i=0;i<3;i++){
                mu_save[ia][i] = mu_ind[ia][i];
            }
        }
        /* -------------------------*/
        induced_moment(mu_ind,ef);
        /* -------------------------*/
        /* correction mu_ind = W mu_ind + (1-W) mu_prev */
        for (int ia=0;ia<nion;ia++){
            for (int i=0;i<3;i++){
                mu_ind[ia][i] = mu_ind[ia][i]*W_ASPC + (1.0 - W_ASPC)*mu_save[ia][i];
            }
        }
    }

    free(mu_save);

}

/******************************************************************************/
/* Root mean square deviation :                                               */
/*          mu = alpha E                                                      */
/*      rmsd = sqrt ((mu/alpha -E)^2/ npol )                                  */
/******************************************************************************/
double get_rmsd_mu(double (*mu)[3], double (*ef)[3]){

    double rmsd=0.0;
    int npol=0;
    int it;

    for (int ia=0;ia<nion;ia++){
        it = typia[ia];
        if (!lpolar[it]) continue;
        npol += 1;
        for (int i=0;i<3;i++){
            if (polia[ia][i][i] != 0.0 ){
                rmsd+= (mu[ia][i] / polia[ia][i][i] - ef[ia][i])*
                       (mu[ia][i] / polia[ia][i][i] - ef[ia][i]);
            }
        }
    }
    return sqrt(rmsd/(double)npol);
}

/******************************************************************************/
/* SCF kO algo                                                                */
/* by Kirill Okhotnikov  somewhere in 2014                                    */
/* the main idea is to perform most of calculation in real space (ewal_dir)   */
/* as the electric dipole are shorter ranged copare to charges                */
/******************************************************************************/
void momentpolaSCFkO(double (*mu_ind)[3],double *upol){

#ifdef DEBUG_SCF
    printf("inside momentpolaSCF \n");
#endif
    int it;
    int iscf;
    double rmsd;
    double (*ef)[3];
    double (*ef_stat)[3];
    double (*ef_ind)[3];
    double (*f_ind)[3];
    double (*dmu_ind)[3];
    double u_coul_stat,u_coul_ind;
    double u_coul_pol;
    double uupol;
    double pdummy;
    double tdummy[3][3];
    double (*efgdummy)[3][3];
    double (*zeroch);
    zeroch=malloc(nion*sizeof(*zeroch));
    ef=malloc(nion*sizeof(*ef));
    efgdummy=malloc(nion*sizeof(*efgdummy));
    ef_stat=malloc(nion*sizeof(*ef_stat));
    ef_ind=malloc(nion*sizeof(*ef_ind));
    f_ind=malloc(nion*sizeof(*f_ind));
    dmu_ind=malloc(nion*sizeof(*dmu_ind));
    for(int ia=0;ia<nion;ia++){
        zeroch[ia]=0.0;
        for(int i=0;i<3;i++){
            ef[ia][i]=0.0;
            ef_stat[ia][i]=0.0;
            ef_ind[ia][i]=0.0;
            f_ind[ia][i]=0.0;
            dmu_ind[ia][i]=0.0;
            for(int j=0;j<3;j++){
                efgdummy[ia][i][j]=0.0;
                tdummy[i][j]=0.0;
            }
        }
    }

    /* static electric field from
     * charges+dipoles+quadrupoles
     * and damping
    */
    multipole_ES(qia,dipia,quadia,&u_coul_stat,&pdummy,tdummy,ef_stat,efgdummy,lqch,ldip,lqua,ldmp,
                /* do   forces stress ef    efg    dir   rec  inpim*/
                        false, false, true, false, true, true, true);

    /* first electric field is static */
    for (int ia=0;ia<nion;ia++){
        for (int i=0;i<3;i++){
            ef[ia][i]=ef_stat[ia][i];
            f_ind[ia][i]=ef_stat[ia][i];
        }
    }
    /* SCF LOOP */
    iscf=0;
    rmsd=DBL_MAX;
    if ( iopnode(istep,npas,nprint) ) {
        printf("  -----------------------------------------------------------\n");
        printf("                 PIM :     self consistent (kOv1)            \n");
        printf("  -----------------------------------------------------------\n");
        printf("        iter          u_pol          u_coul            rmsd  \n");
        printf("  -----------------------------------------------------------\n");
    }

    while ( ( (iscf<max_scf_pim_iter) && (rmsd>conv_tol_ind))  || (iscf<min_scf_pim_iter) ) {

        /* predictor */
        if ( ! iscf ) {
            extrapolate_dipole_aspc(mu_ind,ef,1);
        }
        else {
            induced_moment_inner_kO(dmu_ind,f_ind);
            for (int ia=0;ia<nion;ia++){
                for (int i=0;i<3;i++){
                    mu_ind[ia][i] += dmu_ind[ia][i];
                }
            }
        }
        uupol=get_upol(mu_ind);

        /* induced electric field from only induced dipoles */
        /* quadia -> quad_ind if induce electric field gradient is wanted  */
        /*                                                                   qch   dip  quad  dmp */
        multipole_ES(zeroch,mu_ind,quadia,&u_coul_ind,&pdummy,tdummy,ef_ind,efgdummy,false,true,false,false,
        /* do   forces stress ef    efg    dir   rec  inpim */
                false, false, true, false, true, true, true);

        /* ef = ef_stat + ef_ind */
        for (int ia=0;ia<nion;ia++){
            for (int i=0;i<3;i++){
                ef[ia][i] = ef_stat[ia][i] + ef_ind[ia][i];
            }
        }
        u_coul_pol = u_coul_stat - u_coul_ind;

        /* corrector */
        if (!iscf) extrapolate_dipole_aspc(mu_ind,ef,2);


        rmsd=get_rmsd_mu(mu_ind, ef);
        /* kO */
        for (int ia=0;ia<nion;ia++){
            for (int i=0;i<3;i++){
                it = typia[ia];
                if (!lpolar[it]) continue;
                f_ind[ia][i] = ef[ia][i] - mu_ind[ia][i]*invepolia[ia][i][i];
            }
        }
        /* kO */
        if ( iopnode(istep,npas,nprint ) )  {
            printf("  scf: %5d  "ee3"\n",iscf+1,uupol*coul_unit,u_coul_pol*coul_unit,rmsd);
        }
        iscf+=1;
    }

    //store induced dipole at t
    for (int ia=0;ia<nion;ia++){
        for (int i=0;i<3;i++){
            dipia_ind[ia][0][i]=mu_ind[ia][i];
        }
    }
    if ( iopnode(istep,npas,nprint ) )  {
        putchar('\n');
        printf("  scf converged in %d iterations, rmsd %e\n",iscf,rmsd);
        putchar('\n');
#ifdef DEBUG_ASPC
        printf("  ASPC scf dipia_ind ");
        for (int k=0;k<extrapolate_order;k++){
            printf("%e ",dipia_ind[0][k][0]);
        }
        printf("-> mu_ind %e\n",mu_ind[0][0]);
#endif
    }

    *upol=uupol*coul_unit;
    free(ef);
    free(efgdummy);
    free(ef_stat);
    free(ef_ind);
    free(f_ind);
    free(dmu_ind);
    free(zeroch);

}


/******************************************************************************/
/* this subroutine calculates the induced moment from the total electric
   field and the polarizability tensor
   Basically used in the SCF loop to get induced_moment.

   \f$ \mu_{i,\alpha} =  p_{i,\alpha,\beta} * E_{i,\beta} \f$

   in  : ef electric field vector define at ion position
   out : mu_ind induced electric dipole define at ion position
   polia is the polarizability tensor
 ******************************************************************************/
void induced_moment(double (*mu_ind)[3], double (*ef)[3]){

    int it;
    /*
    ---------------------------------------------------------------
    \mu_{i,\alpha} =  alpha_{i,\alpha,\beta} * E_{i,\beta}
    ---------------------------------------------------------------
    note on units :
    everything are in internal units :
        [ Efield ] = e / A^2
        [ polia  ] = A ^ 3
        [ mu_ind ] = e A
    ---------------------------------------------------------------
    */
    for (int ia=0;ia<nion;ia++){
        it = typia[ia];
        if (!lpolar[it]) continue;
        for (int i=0;i<3;i++){
            mu_ind[ia][i]=0.0;
            for (int j=0;j<3;j++){
                mu_ind[ia][i] += polia[ia][j][i]*ef[ia][j];
            }
        }
    }
}

/******************************************************************************/
/* get induced electric dipole from real part sum only                        */
/* author : kO                                                                */
/******************************************************************************/
void induced_moment_inner_kO(double (*mu_ind)[3], double (*ef_ext)[3]){

#ifdef DEBUG_INNER_KO
    printf("inside induced_moment_inner_kO\n");
#endif
    double udummy;
    int iscf_inner;
    double rmsd_inner,rmsd_ext;
    double alphaES_save;
    iscf_inner=0;
    rmsd_inner=DBL_MAX;
    alphaES_save = alphaES;
    double (*mu_prev)[3],(*zerovec)[3],(*ef_tot)[3],(*ef_ind)[3];
    double (*zeroch);
    double pdummy;
    double tdummy[3][3];
    double (*efgdummy)[3][3];
    zeroch=malloc(nion*sizeof(*zeroch));
    mu_prev=malloc(nion*sizeof(*mu_prev));
    zerovec=malloc(nion*sizeof(*zerovec));
    ef_tot=malloc(nion*sizeof(*ef_tot));
    ef_ind=malloc(nion*sizeof(*ef_ind));
    efgdummy=malloc(nion*sizeof(*efgdummy));
    for(int ia=0;ia<nion;ia++){
        zeroch[ia]=0.0;
        for(int i=0;i<3;i++){
            mu_prev[ia][i]=0.0;
            zerovec[ia][i]=0.0;
            ef_tot[ia][i]=ef_ext[ia][i];
            ef_ind[ia][i]=0.0;
            for(int j=0;j<3;j++){
                efgdummy[ia][i][j]=0.0;
            }

        }
    }
    rmsd_ext = get_rmsd_mu(zerovec, ef_tot);
    induced_moment(mu_ind,ef_tot);

    for(int ia=0;ia<nion;ia++){
        for(int i=0;i<3;i++){
            mu_ind[ia][i]*=omegakO;
        }
    }

    /* inner loop */
    while ((iscf_inner <=10) && (rmsd_inner > 0.1*rmsd_ext) ) {

        iscf_inner+=1;
        /* alpha is reduced for the short-range dipole-dipole interaction
         only real part is calculated in the inner loop*/
        alphaES = 0.001;

        for(int ia=0;ia<nion;ia++){
            for(int i=0;i<3;i++){
                mu_prev[ia][i]=mu_ind[ia][i];
            }
        }

        /* induced electric field from only induced dipoles */
        /* quadia -> quad_ind if induce electric field gradient is wanted  */
        /*                                                                   qch   dip  quad  dmp */
        multipole_ES(zeroch,mu_ind,quadia,&udummy,&pdummy,tdummy,ef_ind,efgdummy,false,true,false,false,
        /* do   forces stress ef    efg    dir   rec   inpim */
                false, false, true, false, true, false, true);

        for(int ia=0;ia<nion;ia++){
            for(int i=0;i<3;i++){
                ef_tot[ia][i] = ef_ext[ia][i] + ef_ind[ia][i];
            }
        }

        rmsd_inner = get_rmsd_mu(mu_ind, ef_tot);
        induced_moment(mu_ind,ef_tot);
        for(int ia=0;ia<nion;ia++){
            for(int i=0;i<3;i++){
                mu_ind[ia][i]= (1-omegakO)*mu_prev[ia][i] + omegakO*mu_ind[ia][i] ;
            }
        }
    }
#ifdef DEBUG_INNER_KO
    if ( iopnode(istep,npas,nprint ) )  {
        printf("  inner loop converged in %d steps   rmsd_inner = %e\n",iscf_inner,rmsd_inner);
    }
#endif
    /* restore alphaES */
    alphaES = alphaES_save;

    free(mu_prev);
    free(zeroch);
    free(zerovec);
    free(ef_tot);
    free(ef_ind);
    free(efgdummy);
#ifdef DEBUG_INNER_KO
    printf("out of induced_moment_inner_kO\n");
#endif
}
