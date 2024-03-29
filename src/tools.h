#ifndef TOOLS_H
#define TOOLS_H
struct sort
{
    double value;
    int index;
};
double dmin(double a,double b);
double dmin3(double a,double b,double c);
double dmin_arr(double *arr, int n);
void expro(double H[3], double U1[3], double U2[3]);
void sliceC(double S[3],double A[3][3],int index);
int sum(int* arr,size_t n);
int check_string(char* label,char buffer[MAX_LEN+1], char allwd[][MAX_LEN+1], int size_allwd);
int check_boolstring(char* label,char buffer[MAX_LEN+1]);
void com(double*,double*,double*,int, double [][3]);
void mainsort(double *arr, int n , int *index);
#endif
