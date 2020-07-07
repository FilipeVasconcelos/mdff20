#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io.h"
#include "nonbonded.h"
/******************************************************************************/
int read_nonbonded(char* controlfn){

    char buffer[MAX_LEN+1];
    FILE * fp;
    fp = fopen (controlfn, "r");
    if (NULL == fp )  {
        pError("opening control file (reading nonbondede)");
        return (-1);
    }
    while (EOF != fscanf(fp, "%s", buffer)) {
        //  lsymmetric
        if (strcmp(buffer,"lsymmetric") == 0 ) {
            lsymmetric=true;
        }
    }
    return 0;
}


void init_nonbonded(char * controlfn){
    read_nonbonded(controlfn);
}
