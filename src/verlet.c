#include "verlet.h"

void allocate_verletlist(VERLETL data){
    data->list=malloc(nion*VNLMAX*sizeof(*data->list));
    data->point=malloc((nion+1)*sizeof(*data->point));
}
void free_verletlist(VERLETL data){
    free(data->list);
    free(data->point);
}

void gen_pbc_verletlist(){
   //verlet_nb.list[ 
}
void check_verletlist(){
}
