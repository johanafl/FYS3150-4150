#include "energy_solver.h"

int main()
{   
    int the_magic_seed = 1572032584;
    int spin_matrix_dim = 2;
    int mc_iterations = 1e7;
    bool convergence = false;
    
    double initial_temp = 1;
    double final_temp = 1;
    double dtemp = 1;       // temperature step length

    time_t seed;
    time(&seed);
    
    IsingModel convergence_model(spin_matrix_dim, mc_iterations, seed);
    convergence_model.iterate_temperature(initial_temp, final_temp, dtemp, convergence);

    return 0;
}