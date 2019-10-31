// #include <mpi.h>
#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include "circular_matrix.cpp"
double const pi = 3.14159265359;

void total_energy_and_magnetization(CircularMatrix& spin, int n,
    double& total_energy, double& total_magnetization)
{   /*
    Computing from upper left in matrix and down to the right.

    Parameters
    ----------
    spin : CircularMatrix
        n x n matrix.
    
    n : int
        Grid dimension.

    total_energy : double refernce
        Total energy.

    total_magnetization : double refernce
        Total magnetic moment.
    */
    for (int i = 0; i < n; i++)
    {   // looping over rows 
        for (int j = 0; j < n; j++)
        {   /*
            looping over columns
            spin(i, j):     spin of the current indices
            spin(i+1, j):   spin below
            spin(i, j+1):   spin to the right
            */
            total_energy        -= spin(i, j, true)*(spin(i, j+1, true) + spin(i+1, j, true));
            total_magnetization += spin(i, j, true);
        }
    }
}

void metropolis_flap(CircularMatrix& spin, double& total_energy,
    double& total_magnetization, int row, int col, double metropolis_random,
    double temperature, double* exp_delta_energy)
{   /*
    Flips a spin and calculates the energy difference. Accepts/rejects the flip
    based on the Metropolis algorithm.

    
    Parameters
    ----------
    spin : CircularMatrix
        Matrix of spin values.
    
    total_energy : double reference
        Total energy of all the spins.

    total_magnetization : double reference
        Total magnetic moment of the system.

    row : int
        A randomly chosen row index.

    col : int
        A randomly chosen column index.

    metropolis_random : double
        Random variable drawn from a uniform distribution on the interval
        [0, 1). Metropolis condition.

    temperature : double
        Temperature of the system.

    exp_delta_energy : double pointer
        Array containing the exponential of the possible energies.
    
    Note
    ----
    spin_here  = -1*spin[row][col]
    spin_left  = spin[row][col-1]
    spin_above = spin[row-1][col]
    spin_right = spin[row][col+1]
    spin_below = spin[row+1][col]
    */
    double spin_here = (-1)*spin(row, col, true);
    double delta_energy = 2*spin_here*(spin(row-1, col, true) + spin(row+1, col, true)
                            + spin(row, col+1, true) + spin(row, col-1, true));

    if (delta_energy <= 0)
    {   // accept new energy if difference is negative        
        spin(row, col, true)*= -1;
        total_energy        += delta_energy;
        total_magnetization += spin_here;
    }
    else if (metropolis_random < exp_delta_energy[(int) (delta_energy + 8)])
    // else if ( (delta_energy > 0) and (metropolis_random < std::exp(delta_energy/temperature)) )
    {   // checks if energy difference is positive and the metropolis condition true
        spin(row, col, true) *= -1;
        total_energy         += delta_energy;
        total_magnetization  += spin_here;
    }
}

int run_shit(int seed)
{   
    // MPI_Init(NULL, NULL);
    // int world_rank;
    // int world_size;
    // MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // MPI_Finalize();

    // int seed = 1337;
    // time_t seed;
    // time(&seed);
    int n = 20;    // grid points
    double temperature = 1;
    double total_magnetization;
    double total_energy;
    int average_runs = 5;           // number of times the spin flip loop is run
    CircularMatrix spin(n, seed);   // initializing matrix with spins
    std::ofstream E_data;
    std::ofstream M_data;

    // std::string filename = "data_files/E_M_data_T=" + std::to_string(temperature) + ".txt";
    E_data.open("data_files/E_data.txt", std::ios_base::app);
    M_data.open("data_files/M_data.txt", std::ios_base::app);
    
    std::mt19937 engine(seed);
    std::uniform_int_distribution<int> uniform_discrete(0, n - 1);
    std::uniform_real_distribution<double> uniform_continuous(0, 1);

    total_energy_and_magnetization(spin, n, total_energy, total_magnetization);

    double J = 1;
    double* exp_delta_energy = new double[17];

    // pre-calculated exponential values
    exp_delta_energy[0]  = std::exp(-8*J/temperature);
    exp_delta_energy[4]  = std::exp(-4*J/temperature);
    exp_delta_energy[8]  = 1;
    exp_delta_energy[12] = std::exp(4*J/temperature);
    exp_delta_energy[16] = std::exp(8*J/temperature);

    for (int i = 0; i < average_runs; i++)
    {   // writing file header
        E_data << std::setw(15) << i;
        M_data << std::setw(15) << i;
    }

    E_data << "\n";
    M_data << "\n";

    for (int j = 0; j < average_runs; j++)
    {   // loops over n*n spin flips a given amount of times
        // saves relevant data for each iteration

        for (int i = 0; i < n*n; i++)
        {   // flips n*n randomly drawn spins in the spin matrix
            
            int row = uniform_discrete(engine);
            int col = uniform_discrete(engine);
            double metropolis_random = uniform_continuous(engine);
            
            metropolis_flap(spin, total_energy, total_magnetization, row, col, metropolis_random, temperature, exp_delta_energy);
        }
        
        E_data << std::setw(15) << temperature;
        E_data << std::setw(15) << total_energy;

        M_data << std::setw(15) << temperature;
        M_data << std::setw(15) << total_magnetization;

    }



    // spin.print();
    std::cout << "tot M: " << total_magnetization << ", tot_E: " << total_energy << std::endl;
    std::cout << seed << std::endl;
    std::cout << std::endl;
    // best value, forgot seed!!
    delete[] exp_delta_energy;

    return 0;
}

int main()
{   
    int the_magic_seed = 1572032584;
    run_shit(the_magic_seed);

    return 0;
}