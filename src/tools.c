#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "io.h"

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

/******************************************************************************/
/* return true if index is even  */
int check_boolstring(char* label,char buffer[MAX_LEN+1]){
    return check_string(label,buffer,allwd_FT_str,16) %2==0;
}
