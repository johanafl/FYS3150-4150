#include "solver.h"


arma::vec acceleration(arma::vec u, double t)
{   
    arma::vec value(3);
    value.zeros();
    double r;
    double x;
    double y;
    double z;

    x = u(0);
    y = u(1);
    z = u(2);
    r = std::sqrt(x*x + y*y + z*z);

    value(0) -= 4*pi*pi*x/(r*r*r);
    value(1) -= 4*pi*pi*y/(r*r*r);
    value(2) -= 4*pi*pi*z/(r*r*r);

    // a_x = -4*pi*pi*pos_x[i]/(r*r*r);
    // a_y = -4*pi*pi*pos_y[i]/(r*r*r);

    return value;
}

arma::vec f(arma::vec u, double t)
{
    return arma::zeros<arma::vec>(2);
}

int main()
{
    int num_steps_input = 1e5;
    double dt_input = 1e-3;

    double init_pos_x = 1;
    double init_pos_y = 0;
    double init_vel_x = 0;
    double init_vel_y = 2*pi;
    
    // system.forward_euler();
    // // system.velocity_verlet();
    // system.write_to_file();

    ForwardEuler lol(num_steps_input);
    lol.set_initial_conditions(init_pos_x, init_pos_y, init_vel_x, init_vel_y);
    lol.solve(acceleration, dt_input);
    lol.write_to_file();


    return 0;
}