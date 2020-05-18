#define VNLMAX 2000

struct VERLETL verlet_nb;
typedef struct VERLETL
{
    int *list,*point; 
    double cut;
    char* label;
}; VERLETL

void allocate_verletlist(VERLETL data);
void free_verletlist(VERLETL data);
void gen_pbc_verletlist();
void check_verletlist();
