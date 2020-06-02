#ifndef INTEGRATION_H
#define INTEGRATION_H
void prop_agate();
void prop_velocity_verlet();
void prop_leap_frog();

void nhcn();
void chain_nhcn(double *kini, double vxi[], double xi[], double Q[], double L);

#endif
