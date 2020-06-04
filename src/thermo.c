#include "stdio.h"
#include "io.h"
#include "constants.h"
#include "thermo.h"
#include "md.h"
#include "config.h"
#include "cell.h"

/******************************************************************************/
void info_thermo(int key, FILE* fp){

    double iso=0.0;
    double l = 3.0 * ((double) nion);
    double temp_r= ( 2.0 *e_kin ) / ( l * boltz_unit );
    double time = ( (double) istep ) * dt / time_unit;
    double acell=simuCell.Anorm[0];
    double bcell=simuCell.Anorm[1];
    double ccell=simuCell.Anorm[2];
    double e_tot=u_lj+e_kin;
    double h_tot = e_tot + e_nvt;
    double pvir_vdw=pvir_lj;
    double pressure=(pvir_lj + temp_r * boltz_unit * simuCell.inveOmega )/press_unit;
    double u_vdw=u_lj;
    double u_tot=u_vdw+u_coul;
    double tau_vdw[3][3];
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            tau_vdw[i][j]=tau_lj[i][j]; /* + other nonbonded contributions */
        }
    }


    if ( key == 0 ) { /* stdout */
        if (ionode) {
            printf("\n");
            printf("  Thermodynamic information                    \n" );
            printf("  ---------------------------------------------\n" );
            printf("  step                  = %9d\n"             ,istep);
            printf("  time                  = "EE"\n"             ,time);
            printf("  Ekin                  = "EE"\n"            ,e_kin);
            printf("  Temp                  = "EE"\n"           ,temp_r);
            printf("  Utot                  = "EE"\n"            ,u_tot);
            printf("  Uvdw                  = "EE"\n"            ,u_vdw);
            printf("  Ucoul                 = "EE"\n"           ,u_coul);
            printf("  Pressure              = "EE"\n"         ,pressure);
            printf("  Pvirvdw               = "EE"\n"         ,pvir_vdw);
            printf("  Pvircoul              = "EE"\n"        ,pvir_coul);
            printf("  volume                = "EE"\n",   simuCell.Omega);
            printf("  a cell                = "EE"\n"            ,acell);
            printf("  b cell                = "EE"\n"            ,bcell);
            printf("  c cell                = "EE"\n"            ,ccell);
            printf("  ---------------------------------------------\n" );
            printf("  Etot                  = "EE"\n"            ,e_tot);
            printf("  Htot                  = "EE"\n"            ,h_tot);
            printf("\n");
            printf("  ---------------------------------------------\n" );
            printf("  non_bonded stress tensor :\n"                    );
            printf("  ---------------------------------------------\n" );
            for (int i=0;i<3;i++){
                printf("  "ee3"\n",tau_vdw[i][0],tau_vdw[i][1],tau_vdw[i][2]);
                iso+=tau_vdw[i][i];
            }
            printf("  ---------------------------------------------\n" );
            printf("  iso = "ee "("ee")\n",iso*onethird,iso*onethirdnion);
            putchar('\n');
        }
    }
    else if (key==1) { /* OSZIFF */
        if (ionode) {
            if (oszcall%10==0) {
                BSEPF(fp);
                fprintf(fp,OSZHEADER);
                BSEPF(fp);
            }
            fprintf(fp,"%7d"ee6"\n",istep,time,e_tot,e_kin,u_tot,u_vdw,u_coul);
            fprintf(fp,"%7d"ee6"\n",istep,temp_r,pressure,pvir_vdw,pvir_coul,simuCell.Omega,h_tot);
            oszcall+=1;
        }
    }
}
