#ifndef PLANET_H
#define PLANET_H

#include "vector_3D.h"

class Planet
{
public:
    double mass;

    double init_x;
    double init_y;
    double init_z;
    Vector_3D pos;
    
    double init_vx;
    double init_vy;
    double init_vz;
    Vector_3D vel;

    Planet();
    Planet(double mass_input, double x, double y, double z, double vx, double vy, double vz);
    double acceleration(double x, double y, double z, double vx, double vy, double vz);
};

#endif