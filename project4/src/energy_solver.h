#ifndef ENERGY_SOLVER_H
#define ENERGY_SOLVER_H

#include "circular_matrix.h"

class IsingModel
{
protected:
// private:
    double* exp_delta_energy = new double[17];  // pre-calculated energy values

    int row;    // row index, will be randomly drawn
    int col;    // column index, will be randomly drawn
    double metropolis_random;   // metropolis condition, will be randomly drawn
    double total_energy;
    double total_magnetization;

    double delta_energy;
    double spin_here;
    
    // values for the averages after convergence
    double sum_total_energy;
    double sum_total_energy_squared;
    double sum_total_magnetization;
    double sum_total_magnetization_absolute;
    double sum_total_magnetization_squared;

    // defining data files
    std::ofstream E_convergence_data;
    std::ofstream M_convergence_data;
    std::ofstream ising_model_data;

    // defining PRNG and distributions
    std::mt19937 engine;
    std::uniform_int_distribution<int> uniform_discrete;
    std::uniform_real_distribution<double> uniform_continuous;


public:
    int n;              // matrix is of dimension nxn
    int mc_iterations;  // number of Monte Carlo iterations
    double J = 1;
    
    // defining the spin matrix
    CircularMatrix spin;

protected:

    void mc_iteration_convergence(double temp);
    void mc_iteration_stable(double temp);
    void iterate_spin_flip(double temp);

public:

    IsingModel(int spin_mat_dim, int mc_iterations_input, long seed);
    void iterate_temperature(double initial_temp, double final_temp,
        double dtemp, bool convergence);
    void total_energy_and_magnetization(CircularMatrix& spin, int n,
        double& total_energy, double& total_magnetization);
    void set_new_input(int spin_mat_dim, int mc_iterations_input, double inter_strenght_J, long seed);
    void set_interactions_strength(double strength_J);
    void set_mc_iterations(int mc_iterations_input);
    void set_spin_dim(int spin_mat_dim);
    void set_order_spins();
    void iterate_monte_carlo_cycles(int initial_MC, int final_MC, int dMC);
    void metropolis_flap(CircularMatrix& spin, double& total_energy,
        double& total_magnetization, int row, int col, double metropolis_random,
        double temperature, double* exp_delta_energy);
    ~IsingModel();
};


#endif