#include <stdio.h> 
#include <stdlib.h> 
#include <time.h> 
#include <math.h> 
  
#include "constants.h"

void init_rand(time_t t)
{
    srand((unsigned) time(&t));
}

// Generates random numbers in range [lower, upper]. 
double randD() 
{
    return (double) rand()/(double)(RAND_MAX); 
    
} 
  

double box_muller(double mean, double sigma)
{
    double C,U,V,R,G;
    U=randD();
    V=randD();
    R = sqrt( -2.0 * log(U));
    C = cos ( TPI * V );
    G = R * C;
    return mean + G * sigma;

}
