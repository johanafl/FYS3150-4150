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
{   /*
    One step with Forward Euler.

    Parameters
    ----------
    f : function pointer
        The right hand side of the ODE.

    k : int
        Current step in the integration.
    */

    arma::vec a = f(pos.col(k), 0);

    pos.col(k+1) = pos.col(k) + dt*vel.col(k);
    vel.col(k+1) = vel.col(k) + dt*a;    
}

void VelocityVerlet::advance(arma::vec (*f)(arma::vec u, double t), int k)
{   /*
    One step with Velocity Verlet.

    Parameters
    ----------
    f : function pointer
        The right hand side of the ODE.

    k : int
        Current step in the integration.
    */
    
    arma::vec a1 = f(pos.col(k), 0);
    pos.col(k+1) = pos.col(k) + dt*vel.col(k) + dt*dt/2*a1;

    arma::vec a2 = f(pos.col(k+1), 0);
    vel.col(k+1) = vel.col(k) + dt/2*(a2 + a1);

}