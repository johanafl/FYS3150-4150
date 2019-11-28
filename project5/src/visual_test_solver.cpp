#include "solver.h"

int main()
{
    int num_steps_input = 1e5;
    double dt_input = 1e-3;

    double init_pos_x = 1;
    double init_pos_y = 0;
    double init_vel_x = 0;
    double init_vel_y = 2*pi;
    
    Solver system(num_steps_input, dt_input, init_pos_x, init_pos_y, init_vel_x, init_vel_y);
    // system.forward_euler();
    system.velocity_verlet();
    system.write_to_file();
    
    return 0;
}