#ifndef RDF_H
#define RDF_H
/* ****************************************************************************/
/*                              global parameters                             */
/* ****************************************************************************/
double                               resg; /* resolution */
int                                 nconf; /* number of configuration */
double                              cutgr; /* */
int                                 nbins; /* */
int            (*gr)[NTYPEMAX][NTYPEMAX]; /* */
/* ****************************************************************************/
/*                                prototypes                                  */
/* ****************************************************************************/
int read_rdf(char* controlfn);
void init_rdf(char* controlfn);
void check_rdf();
void info_rdf();
void run_rdf();
void alloc_rdf();
void free_rdf();
void get_gr();
int  write_rdf();
int rdf_readtraj();
#endif
