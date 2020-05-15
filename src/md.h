#ifndef MD_H
#define MD_H
double          temp;
int             npas;
int           nprint;
int           nequil;
double            dt;
double tauTberendsen;
int            itime;

int read_md(char* controlfn);
void init_md(char* controlfn);
void info_md();
void run_md();
#endif
