#include "stdio.h"
#include "io.h"
#include "constants.h"
#include "thermo.h"
#include "md.h"
#include "config.h"
#include "cell.h"

void write_thermo(){

    double time = ( (double) itime ) * dt / time_unit;
    double omega=simu_cell.Omega;
    double acell=simu_cell.Anorm[0];
    double bcell=simu_cell.Anorm[1];
    double ccell=simu_cell.Anorm[2];
    double e_tot=u_lj+e_kin;

    if (ionode) {
        printf("\n");
        printf("  Thermodynamic information                    \n");
        printf("  ---------------------------------------------\n");
        printf("  step                  = %9d\n"            ,itime);
        printf("  time                  = %19.12e\n"             ,time);
        printf("  Ekin                  = %19.12e\n"        ,e_kin);
        printf("  Temp                  = %19.12e\n"       ,temp_r);
        printf("  Utot                  = %19.12e\n"         ,u_lj);
        printf("  Pressure              = xxxxxxxxxx           \n");
        printf("  volume                = %19.12e\n"        ,omega);
        printf("  a cell                = %19.12e\n"        ,acell);
        printf("  b cell                = %19.12e\n"        ,bcell);
        printf("  c cell                = %19.12e\n"        ,ccell);
        printf("  ---------------------------------------------\n");
        printf("  Etot                  = %19.12e\n"        ,e_tot);
        printf("  Htot                  = %19.12e\n"        ,e_tot);
        printf("\n");
        printf("\n");
        printf("non_bonded stress tensor\n");
        double iso=0.0;
        for (int i=0;i<3;i++){
            printf("%15.8e %15.8e %15.8e\n",tau_nonb[i][0],tau_nonb[i][1],tau_nonb[i][2]);
            iso+=tau_nonb[i][i];
        }
        printf("iso = %15.8e (%15.8e)\n",iso*onethird,iso*onethirdnion);
        putchar('\n');

    }
}
