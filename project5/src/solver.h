#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
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

    Solver(int num_steps_input, double dt_input, double init_pos_x, 
        double init_pos_y, double init_vel_x, double init_vel_y);
    void one_step_forward_euler(double dt, double pos_x, double pos_y, 
        double vel_x, double vel_y, double a_x, double a_y, double& pos_x_new, 
        double& pos_y_new, double& vel_x_new, double& vel_y_new);
    void forward_euler();
    void one_step_velocity_verlet(double dt, double abs_dist, double pos_x, 
        double pos_y, double vel_x, double vel_y, double a_x, double a_y, 
        double& pos_x_new, double& pos_y_new, double& vel_x_new, 
        double& vel_y_new);
    void velocity_verlet();
    void write_to_file();
    ~Solver();
};

#endif