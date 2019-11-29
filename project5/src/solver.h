#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <armadillo>

const double pi = 3.14159265358979323846;

template <class T>
class Solver
{
public:
    int num_steps;
    int num_stellar_objects;
    
    double dt;

    double* pos_x;
    double* pos_y;
    double* vel_x;
    double* vel_y;

    arma::mat vel;
    arma::mat pos;

    Solver(int num_steps_input, int num_stellar_objects);

    // Virtual to allow a method in the superclass to call a method advance
    // in a subclass.
    virtual void advance(T object, int k);
    
    void set_initial_conditions(arma::vec U0);
    
    void solve(T object, double dt_input);
    
    void write_to_file();
    
    ~Solver();
};

template <class T>
class ForwardEuler : public Solver<T>
{
public:
    using Solver<T>::Solver;   // For inheriting the constructor of Solver.
    // void advance(arma::vec (*f)(arma::vec u, double t), int k);
    void advance(T object, int k); 
};

template <class T>
class VelocityVerlet : public Solver<T>
{
public:
    using Solver<T>::Solver;   // For inheriting the constructor of Solver.
    void advance(T object, int k);
};

#endif