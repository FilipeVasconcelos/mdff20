#ifndef VERLET_H
#define VERLET_H
#define VNLMAX 500


double *xvl,*yvl,*zvl; /* last positions when verlet list was updated */  
int updatevl;

struct VERLETL *verlet_nb;

typedef struct VERLETL
{
    int *list,*point; 
    double cut;
    char* label;
} VERLETL;

VERLETL* allocate_verletlist(char* label);
void free_verletlist(char* label);
void gen_pbc_verletlist();
void check_verletlist();
#endif
