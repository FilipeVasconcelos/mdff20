// Structure atom decomposition parallelism
typedef struct DEC
{
    int dimData;
    int iaStart;
    int iaEnd;
    char* label;
}DEC ;

void do_split(int n,int np,int mrank,DEC dec,char* lab);
