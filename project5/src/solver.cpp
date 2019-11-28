#include "solver.h"

Solver::Solver(int num_steps_input)
    : pos(3, num_steps_input+1), vel(3, num_steps_input+1)
{
    num_steps = num_steps_input;
}

void Solver::set_initial_conditions(double init_pos_x, double init_pos_y,
        double init_vel_x, double init_vel_y)
{
    // pos_x[0] = init_pos_x;
    // pos_y[0] = init_pos_y;
    // vel_x[0] = init_vel_x;
    // vel_y[0] = init_vel_y;
    double init_pos_z = 0;
    double init_vel_z = 0;

    pos(0, 0) = init_pos_x;
    pos(1, 0) = init_pos_y;
    pos(2, 0) = init_pos_z;
    vel(0, 0) = init_vel_x;
    vel(1, 0) = init_vel_y;
    vel(2, 0) = init_vel_z;
    
}

void Solver::solve(arma::vec (*f)(arma::vec u, double t), double dt_input)
{
    dt = dt_input;
    double t = 0;  // Currently unused.
    
    for (int k = 0; k < num_steps; k++)
    {
        advance(f, k);
    }
}

void Solver::advance(arma::vec (*f)(arma::vec u, double t), int k)
{
    // Dummy method.
    std::cout << "NotImplementedError" << std::endl;
}


void Solver::one_step_velocity_verlet(double dt, double abs_dist, double pos_x, 
    double pos_y, double vel_x, double vel_y, double a_x, double a_y, 
    double& pos_x_new, double& pos_y_new, double& vel_x_new, 
    double& vel_y_new)
{
    // Velocity verlet
    pos_x_new = pos_x + dt*vel_x + dt*dt/2*a_x;
    pos_y_new = pos_y + dt*vel_y + dt*dt/2*a_y;
    // pos_x[i+1] = pos_x[i] + dt*vel_x[i] + dt*dt/2*a_x;
    // pos_y[i+1] = pos_y[i] + dt*vel_y[i] + dt*dt/2*a_y;

    double a_x_next = -4*pi*pi*pos_x_new/(abs_dist*abs_dist*abs_dist);
    double a_y_next = -4*pi*pi*pos_y_new/(abs_dist*abs_dist*abs_dist);
    // a_x_next = -4*pi*pi*pos_x[i+1]/(r*r*r);
    // a_y_next = -4*pi*pi*pos_y[i+1]/(r*r*r);
    
    vel_x_new = vel_x + dt/2*(a_x_next + a_x);
    vel_y_new = vel_y + dt/2*(a_y_next + a_y);
    // vel_x[i+1] = vel_x[i] + dt/2*(a_x_next + a_x);
    // vel_y[i+1] = vel_y[i] + dt/2*(a_y_next + a_y);            
}

void Solver::velocity_verlet()
{
    double a_x;
    double a_y;
    double a_x_next;
    double a_y_next;
    double r;

    // Verlet
    for (int i = 0; i<num_steps-1; i++)
    {
        // std::cout << std::endl;
        r = sqrt(pos_x[i]*pos_x[i] + pos_y[i]*pos_y[i]);

        a_x = -4*pi*pi*pos_x[i]/(r*r*r);
        a_y = -4*pi*pi*pos_y[i]/(r*r*r);

        one_step_velocity_verlet(dt, r, pos_x[i], pos_y[i], vel_x[i], 
            vel_y[i], a_x, a_y, pos_x[i+1], pos_y[i+1], vel_x[i+1], vel_y[i+1]);

        // pos_x[i+1] = pos_x[i] + dt*vel_x[i] + dt*dt/2*a_x;
        // pos_y[i+1] = pos_y[i] + dt*vel_y[i] + dt*dt/2*a_y;

        // a_x_next = -4*pi*pi*pos_x[i+1]/(r*r*r);
        // a_y_next = -4*pi*pi*pos_y[i+1]/(r*r*r);
        
        // vel_x[i+1] = vel_x[i] + dt/2*(a_x_next + a_x);
        // vel_y[i+1] = vel_y[i] + dt/2*(a_y_next + a_y);
    }
}

void Solver::write_to_file()
{
    // defining data files
    std::ofstream tull_mc_tull;
    // tull_mc_tull.open("data_files/tull_mc_tull.txt", std::ios_base::app);
    tull_mc_tull.open("data_files/tull_mc_tull.txt");
    
    for (int i = 0; i < num_steps; i++)
    {
        tull_mc_tull << std::setw(20) << std::setprecision(15);
        tull_mc_tull << dt*i;
        tull_mc_tull << std::setw(30) << std::setprecision(15);
        tull_mc_tull << pos(0, i);
        tull_mc_tull << std::setw(30) << std::setprecision(15);
        tull_mc_tull << pos(1, i);
        tull_mc_tull << std::setw(30) << std::setprecision(15);
        tull_mc_tull << vel(0, i);
        tull_mc_tull << std::setw(30) << std::setprecision(15);
        tull_mc_tull << vel(1, i);
        tull_mc_tull << std::endl;
    }
    tull_mc_tull.close();
}

Solver::~Solver()
{
}


void ForwardEuler::advance(arma::vec (*f)(arma::vec u, double t), int k)
{   
    arma::vec a = f(pos.col(k), 0);

    vel(0, k+1) = vel(0, k) + dt*a(0);
    vel(1, k+1) = vel(1, k) + dt*a(1);

    pos(0, k+1) = pos(0, k) + dt*vel(0, k);
    pos(1, k+1) = pos(1, k) + dt*vel(1, k);
    
}