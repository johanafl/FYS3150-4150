#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include "solver.h"

class Solarsystem
{
private:
    
    int num_planets = 0;
    int array_size  = 3;
    double* pos  = new double[array_size];
    double* vel  = new double[array_size];
    double* mass = new double[array_size/3];

public:
    Solarsystem();    
    void add_planet(double mass_input, arma::vec U0);
    arma::vec get_U0();
    void resize();
    arma::vec acceleration(arma::vec u, double t);
    void solve_system();
};

#endif