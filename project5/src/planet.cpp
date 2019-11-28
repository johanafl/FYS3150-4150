#include "planet.h"

Planet::Planet() : pos(0, 0, 0), vel(0, 0, 0){}

Planet::Planet(double mass_input, double x, double y, double z, double vx, 
    double vy, double vz)
    : pos(x, y, z), vel(vx, vy, vz)
{
    mass    = mass_input;

    init_x  = x;
    init_y  = y;
    init_z  = z;

    init_vx = vx;
    init_vy = vy;
    init_vz = vz;
}

double Planet::acceleration(double x, double y, double z, double vx, double vy, double vz)
{
    // return force or acceleration.
}