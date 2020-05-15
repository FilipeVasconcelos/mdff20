#include<stdio.h>

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

void sliceC(double S[3], double A[3][3],int index)
{
    for (int i=0;i<3;i++)
    {
        S[i]=A[index][i];
    }
}

int sum(int *arr,size_t n){
    int s=0; size_t i;
    for(i=0;i<n;i++) s+=arr[i];
    return s;
}
/*void printV3(double V[3])
{
}
*/
