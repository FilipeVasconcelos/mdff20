#ifndef KINETIC_H
#define KINETIC_H
void init_velocities();
void print_velocities();
void maxwellboltzmann_velocities();
void calc_temp(double *temp, double *ekin,int flag);
void rescale_velocities();
#endif
