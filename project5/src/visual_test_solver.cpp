#include "solver.h"
// #include "solarsystem.cpp"


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

    return value;
}

arma::vec f(arma::vec u, double t)
{
    return arma::zeros<arma::vec>(2);
}

int main()
{
    int num_steps_input = 1e5;
    int num_stellar_objects_input = 1;
    double dt = 1e-3;
     

    double init_pos_x = 1;
    double init_pos_y = 0;
    double init_pos_z = 0;
    double init_vel_x = 0;
    double init_vel_y = 2*pi;
    double init_vel_z = 0;

    arma::vec U0 = {init_pos_x, init_pos_y, init_pos_z,
        init_vel_x, init_vel_y, init_vel_z};

    // Solarsystem system();
    // system.add_planet();

    // ForwardEuler solved(num_steps_input, num_stellar_objects_input);
    // // VelocityVerlet solved(num_steps_input);
    // solved.set_initial_conditions(U0);
    // solved.solve(acceleration, dt);
    // solved.write_to_file();

    


    return 0;
}