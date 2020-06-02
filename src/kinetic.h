#ifndef KINETIC_H
#define KINETIC_H
void init_velocities();
void print_velocities();
void maxwellboltzmann_velocities();
double calc_kin();
double calc_temp(double kin);
void rescale_velocities();
#endif
