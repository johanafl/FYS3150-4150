#include "solver.h"

Solver::Solver(int num_steps_input, double dt_input, double init_pos_x, 
    double init_pos_y, double init_vel_x, double init_vel_y)
{
    num_steps = num_steps_input;
    dt = dt_input;

    pos_x = new double[num_steps];
    pos_y = new double[num_steps];
    vel_x = new double[num_steps];
    vel_y = new double[num_steps];

    pos_x[0] = init_pos_x;
    pos_y[0] = init_pos_y;
    vel_x[0] = init_vel_x;
    vel_y[0] = init_vel_y;
}

void Solver::one_step_forward_euler(double dt, double pos_x, double pos_y, 
    double vel_x, double vel_y, double a_x, double a_y, double& pos_x_new, 
    double& pos_y_new, double& vel_x_new, double& vel_y_new)
{
    // Forward Euler
    vel_x_new = vel_x + dt*a_x;
    vel_y_new = vel_y + dt*a_y;
    pos_x_new = pos_x + dt*vel_x;
    pos_y_new = pos_y + dt*vel_y;
}

//void func ( void (*f)(int) ); // shold let the function f pass as a parameter.
void Solver::forward_euler()
{
    double a_x;
    double a_y;
    double r;

    for (int i = 0; i<num_steps-1; i++)
    {
        r = sqrt(pos_x[i]*pos_x[i] + pos_y[i]*pos_y[i]);

        a_x = -4*pi*pi*pos_x[i]/(r*r*r);
        a_y = -4*pi*pi*pos_y[i]/(r*r*r);

        one_step_forward_euler(dt, pos_x[i], pos_y[i], vel_x[i], vel_y[i], 
            a_x, a_y, pos_x[i+1], pos_y[i+1], vel_x[i+1], vel_y[i+1]);

        // pos_x[i+1] = pos_x[i] + dt*vel_x[i];
        // pos_y[i+1] = pos_y[i] + dt*vel_y[i];

        // vel_x[i+1] = vel_x[i] + dt*a_x;
        // vel_y[i+1] = vel_y[i] + dt*a_y;
    }
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
    
    for (int i=0; i<num_steps; i++)
    {
        tull_mc_tull << std::setw(20) << std::setprecision(15);
        tull_mc_tull << dt*i;
        tull_mc_tull << std::setw(30) << std::setprecision(15);
        tull_mc_tull << pos_x[i];
        tull_mc_tull << std::setw(30) << std::setprecision(15);
        tull_mc_tull << pos_y[i];
        tull_mc_tull << std::setw(30) << std::setprecision(15);
        tull_mc_tull << vel_x[i];
        tull_mc_tull << std::setw(30) << std::setprecision(15);
        tull_mc_tull << vel_y[i];
        tull_mc_tull << std::endl;
    }
    tull_mc_tull.close();
}

Solver::~Solver()
{
    delete[] pos_x;
    delete[] pos_y;
    delete[] vel_x;
    delete[] vel_y;
}