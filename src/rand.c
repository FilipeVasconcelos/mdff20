#include <stdio.h> 
#include <stdlib.h> 
#include <time.h> 
#include <math.h> 
  
#include "constants.h"

/******************************************************************************/
void init_rand(time_t t) {
    srand((unsigned) time(&t));
}

/******************************************************************************/
double randD() {
    return (double) rand()/(double)(RAND_MAX); 
} 
  

/******************************************************************************/
double box_muller(double mean, double sigma) {
    
    double C,U,V,R;
    U=randD();
    V=randD();
    while ( (U==0.0) || (V==0.0)) {
        U=randD();
        V=randD();
    }
    R = sqrt( -2.0 * log(U));
    C = cos ( 2.0 * PI * V );
    return mean + R * C * sigma;
}
