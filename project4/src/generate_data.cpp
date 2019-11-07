#include "energy_solver.h"
#include <chrono>

int main()
{   

    // average_time += comp_time.count();
    int the_magic_seed = 1572032584;
    int spin_matrix_dim = 10;
    int mc_iterations = 1e5;
    bool convergence = false;
    
    double initial_temp = 1;
    double final_temp = 1;
    double dtemp = 1;       // temperature step length
    
    int initial_MC = 1e4;
    int final_MC   = 1e5;
    int dMC        = 1e4;

    time_t seed;
    time(&seed);
    
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    
    IsingModel convergence_model(spin_matrix_dim, mc_iterations, seed);
    // convergence_model.iterate_temperature(initial_temp, final_temp, dtemp, convergence);
    convergence_model.iterate_monte_carlo_cycles(initial_MC, final_MC, dMC);
    
    // ending timer
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> comp_time  = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);

    // std::cout << "iterations: " << mc_iterations;
    std::cout << " time: " << comp_time.count() << std::endl;

    return 0;
}