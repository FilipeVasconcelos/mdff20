#include "stdio.h"
#include "io.h"
#include "constants.h"
#include "thermo.h"
#include "md.h"
#include "config.h"
#include "cell.h"

void info_thermo(){

    double iso=0.0;
    double time = ( (double) istep ) * dt / time_unit;
    double acell=simuCell.Anorm[0];
    double bcell=simuCell.Anorm[1];
    double ccell=simuCell.Anorm[2];
    double e_tot=u_lj+e_kin;
    double pvir_lj= vir_lj * simuCell.inveOmega;
    double pressure=(pvir_lj + temp_r * boltz_unit * simuCell.inveOmega )/press_unit;

    if (ionode) {
        printf("\n");
        printf("  Thermodynamic information                    \n" );
        printf("  ---------------------------------------------\n" );
        printf("  step                  = %9d\n"             ,istep);
        printf("  time                  = %19.12e\n"          ,time);
        printf("  Ekin                  = %19.12e\n"         ,e_kin);
        printf("  Temp                  = %19.12e\n"        ,temp_r);
        printf("  U_lj                  = %19.12e\n"          ,u_lj);
        printf("  Pressure              = %19.12e\n"      ,pressure);
        printf("  Pvir_lj               = %19.12e\n"       ,pvir_lj);
        printf("  volume                = %19.12e\n",simuCell.Omega);
        printf("  a cell                = %19.12e\n"         ,acell);
        printf("  b cell                = %19.12e\n"         ,bcell);
        printf("  c cell                = %19.12e\n"         ,ccell);
        printf("  ---------------------------------------------\n" );
        printf("  Etot                  = %19.12e\n"         ,e_tot);
        printf("  Htot                  = %19.12e\n"         ,e_tot);
        printf("\n");
        printf("  ---------------------------------------------\n" );
        printf("  non_bonded stress tensor :\n"                    );
        printf("  ---------------------------------------------\n" );
        for (int i=0;i<3;i++){
            printf("  %14.7e %14.7e %14.7e\n",tau_nonb[i][0],tau_nonb[i][1],tau_nonb[i][2]);
            iso+=tau_nonb[i][i];
        }
        printf("  ---------------------------------------------\n" );
        printf("  iso = %15.7e (%15.7e)\n",iso*onethird,iso*onethirdnion);
        putchar('\n');
    }
}
