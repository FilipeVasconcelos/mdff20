// Structure atom decomposition parallelism
struct DEC atom_dec, * patom_dec;

typedef struct DEC
{
    int dimData;
    int iaStart;
    int iaEnd;
    char* label;
}DEC ;

int do_split(int n,int np,int mrank,DEC* dec,char* lab);
