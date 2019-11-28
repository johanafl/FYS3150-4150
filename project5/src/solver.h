#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <armadillo>
const double pi = 3.14159265358979323846;


class Solver
{
public:
    int num_steps;
    double dt;

    double* pos_x;
    double* pos_y;
    double* vel_x;
    double* vel_y;

    arma::mat vel;
    arma::mat pos;

    Solver(int num_steps_input);

    virtual void advance(arma::vec (*f)(arma::vec u, double t), int k);
    
    void set_initial_conditions(double init_pos_x, double init_pos_y,
        double init_vel_x, double init_vel_y);
    
    void one_step_forward_euler(double dt, double pos_x, double pos_y, 
        double vel_x, double vel_y, double a_x, double a_y, double& pos_x_new, 
        double& pos_y_new, double& vel_x_new, double& vel_y_new);
    
    void solve(arma::vec (*f)(arma::vec u, double t), double dt_input);
    
    void forward_euler();
    
    void one_step_velocity_verlet(double dt, double abs_dist, double pos_x, 
        double pos_y, double vel_x, double vel_y, double a_x, double a_y, 
        double& pos_x_new, double& pos_y_new, double& vel_x_new, 
        double& vel_y_new);
    
    void velocity_verlet();
    
    void write_to_file();
    
    ~Solver();
};

class ForwardEuler : public Solver
{
public:
    using Solver::Solver;
    void advance(arma::vec (*f)(arma::vec u, double t), int k);
};

#endif