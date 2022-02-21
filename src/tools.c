#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "io.h"
#include "config.h"
#include "float.h"
#include "tools.h"

/*******************************************************************************
 EXPRO
 calculates the x-product of two vectors
 adapted from VASP ;)
 ******************************************************************************/
void expro(double H[3], double * U1, double * U2)
{
  H[0]=U1[1]*U2[2]-U1[2]*U2[1];
  H[1]=U1[2]*U2[0]-U1[0]*U2[2];
  H[2]=U1[0]*U2[1]-U1[1]*U2[0];
}

/******************************************************************************/
void sliceC(double S[3], double A[3][3],int index)
{
    for (int i=0;i<3;i++)
    {
        S[i]=A[index][i];
    }
}

/******************************************************************************/
int sum(int *arr,size_t n){
    int s=0; size_t i;
    for(i=0;i<n;i++) s+=arr[i];
    return s;
}
/******************************************************************************/
double dmin(double a,double b){
    return a<b?a:b;
}
/******************************************************************************/
double dmin3(double a,double b,double c){
    return dmin(dmin(a,b),c);
}

/******************************************************************************/
double dmin_arr(double *arr, int n){
    double min=DBL_MAX;
    for (int i=0;i<n;i++){
        if ( arr[i]<min ) min=arr[i];
    }
    return min;
}

/******************************************************************************/
int check_string(char* label,char buffer[MAX_LEN+1], char allwd[][MAX_LEN+1], int size_allwd){
    int check=-1;
    for (int k=0;k<size_allwd;k++){
        if (strcmp(buffer,allwd[k]) == 0){
            check=k;
        }
    }
    if (check<0){
        pError("value not allowed in input file\n");
        printf("%s %s\n",label,buffer);
        exit(-1);
    }
    return check;
}

//int check_fscanf(char* label, int out, int check){
//    return out = check ;
//}


/******************************************************************************/
/* return true if index is even  */
int check_boolstring(char* label,char buffer[MAX_LEN+1]){
    return check_string(label,buffer,allwd_FT_str,16) %2==0;
}

/******************************************************************************/
/* center of mass by type (e.g comit)*/
void com(double *ax,double *ay,double *az, int n, double comit[n][3]){

    for(int it=0;it<ntype+1;it++){
        for(int k=0;k<3;k++){
            comit[it][k]=0.0;
        }
    }
    int it;
    for (int ia=0;ia<nion;ia++){
        it=typia[ia];
        comit[it][0]+=ax[ia];
        comit[it][1]+=ay[ia];
        comit[it][2]+=az[ia];
        comit[ntype][0]+=ax[ia];
        comit[ntype][1]+=ay[ia];
        comit[ntype][2]+=az[ia];
    //    printf("ia %d %15.8e %15.8e %15.8e %15.8e\n",ia,ax[ia],ay[ia],az[ia],comit[ntype][2]);
    }
    for(int it=0;it<ntype+1;it++){
        for(int k=0;k<3;k++){
            comit[it][k]*=invenionit[it];
        //    printf("it %d k %d inve %f %15.8e\n",it,k,invenionit[it],comit[it][k]);
        }
    }
}

/******************************************************************************/
int cmp(const void *a,const void *b){
    struct sort *a1 = (struct sort *)a;
    struct sort *a2 = (struct sort *)b;
    if ((*a1).value < (*a2).value){
        return -1;
    }
    else if ((*a1).value > (*a2).value){
        return 1;
    }
    else {
        return 0;
    }
}

void mainsort(double *arr, int n , int *index) {
    struct sort obj[n];
    for (int i=0;i<n;i++){
        obj[i].value = arr[i];
        obj[i].index = i;
    }
    qsort(obj,n,sizeof(obj[0]),cmp);
    for (int i=0;i<n;i++){
        index[i]=obj[i].index;
    }
}



/******************************************************************************/
void merge(double *a,int n,int m) {
    int i,j,k;
    int *x =malloc(n*sizeof(double));
    for (i=0,j=m,k=0;k<n;k++){
        x[k] = j == n      ? a[i++]
             : i == m      ? a[j++]
             : a[j] < a[i] ? a[j++]
             :               a[i++];
    }
    for (i=0;i<n;i++){
        a[i]=x[i];
    }
    free(x);
}
/******************************************************************************/
void merge_sort(double *a,int n){
    if (n<2) return;
    int m = n/2;
    merge_sort(a,m);
    merge_sort(a+m,n-m);
    merge_sort(a,m);
    merge(a,n,m);
}
