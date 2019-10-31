// #include <mpi.h>
#include "circular_matrix.h"

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
    {   // checks if energy difference is positive and the metropolis condition true
        spin(row, col, true) *= -1;
        total_energy         += delta_energy;
        total_magnetization  += spin_here;
    }
}


class IsingModel
{

private:
    int n;    // grid points
    int MC_iterations = 5;   // number of times the spin flip loop is run
    double* exp_delta_energy = new double[17];
    double J = 1;

    int row;
    int col;
    double metropolis_random;
    double total_energy;
    double total_magnetization;
    
    // setting temperature data
    double initial_temp = 1;
    double final_temp = 10;
    double dtemp = 1;

    int seed;

    std::ofstream E_data;
    std::ofstream M_data;

    // std::mt19937 engine(int seed);
    // std::uniform_int_distribution<int> uniform_discrete(0, n - 1);
    // std::uniform_real_distribution<double> uniform_continuous(0, 1);
    std::mt19937 engine;
    std::uniform_int_distribution<int> uniform_discrete;
    std::uniform_real_distribution<double> uniform_continuous;
    
    CircularMatrix spin;

public:
    IsingModel(int n)
    {
        // seed = 1337;
        // engine(seed);
        // uniform_discrete(0, n-1);
        // uniform_continuous(0, 1);

        spin(4);   // initializing matrix with spins
        
        total_energy_and_magnetization(spin, n, total_energy, total_magnetization);

        // initialising data files
        E_data.open("data_files/E_data.txt", std::ios_base::app);
        M_data.open("data_files/M_data.txt", std::ios_base::app);

        E_data << std::setw(15) << " ";
        M_data << std::setw(15) << " ";

        for (int i = 1; i <= MC_iterations; i++)
        {   // writing file header
            E_data << std::setw(15) << i*n*n;
            M_data << std::setw(15) << i*n*n;
        }

        E_data << "\n";
        M_data << "\n";
    }


    void MC_iteration()
    {   /*
        Runs the spin flip a given amount of times. Generates data for finding
        how many iterations is needed for convergence. Keeps the energy values
        without averaging. Only runs for a single temperature value.
        */

        double single_temp = 1;    // runs the convergence check for a single temperature

        exp_delta_energy[0]  = std::exp(-8*J/single_temp);
        exp_delta_energy[4]  = std::exp(-4*J/single_temp);
        exp_delta_energy[8]  = 1;
        exp_delta_energy[12] = std::exp(4*J/single_temp);
        exp_delta_energy[16] = std::exp(8*J/single_temp);

        for (int j = 0; j < MC_iterations; j++)
        {   // loops over n*n spin flips a given amount of times
            // saves relevant data for each iteration

            iterate_spin_flip(single_temp);
            // writing calculated data to file  
            E_data << std::setw(15) << total_energy;
            M_data << std::setw(15) << total_magnetization;

        }
    }

    void MC_iteration_average(double temp)
    {   /*
        Runs the spin flip a given amount of times. Does not keep every energy
        value, but calculates the average.

        Parameters
        ----------
        temp : double
            Temperature value.
        */

        for (int j = 0; j < MC_iterations; j++)
        {   // loops over n*n spin flips a given amount of times

            // add array business to calculate average
            iterate_spin_flip(temp);

        }

    }

    void iterate_spin_flip(double temp)
    {   /*
        Picks a random row and a random column. Picks a random number for the
        Metropolis condition. Flips the spin with the function metropolis_flap.
        This is done n*n times.

        Parameters
        ----------
        temp : double
            Temperature value for which to calculate the data.
        */
        
        for (int i = 0; i < n*n; i++)
        {   // flips n*n randomly drawn spins in the spin matrix
            
            row = uniform_discrete(engine);
            col = uniform_discrete(engine);
            metropolis_random = uniform_continuous(engine);
            
            metropolis_flap(spin, total_energy, total_magnetization, row, col, metropolis_random, temp, exp_delta_energy);
        }
    }

    void iterate_temperature()
    {   /*
        Iterates over a given set of temperature values.
        */
    
        for (double temp = initial_temp; temp <= final_temp; temp += dtemp)
        {   // looping over temperature values

            // pre-calculated exponential values
            exp_delta_energy[0]  = std::exp(-8*J/temp);
            exp_delta_energy[4]  = std::exp(-4*J/temp);
            exp_delta_energy[8]  = 1;
            exp_delta_energy[12] = std::exp(4*J/temp);
            exp_delta_energy[16] = std::exp(8*J/temp);
            
            // writing temperature values in the first column
            E_data << std::setw(15) << temp;
            M_data << std::setw(15) << temp;
            
            E_data << "\n";
            M_data << "\n";
        }
    }

    ~IsingModel()
    {
        delete[] exp_delta_energy;
    }

};


int generate_data(int seed)
{   /*
    Generates energy and magnetization data for a given set of temperature
    values and a given set of average runs. Writes data to file.

    Parameters
    ----------
    seed : int
        Seed for the Mersenne Twister 19937 PRNG.
    */

    // MPI_Init(NULL, NULL);
    // int world_rank;
    // int world_size;
    // MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // MPI_Finalize();

    // int seed = 1337;
    // time_t seed;
    // time(&seed);

    // std::string filename = "data_files/E_M_data_T=" + std::to_string(temp) + ".txt";








    // // spin.print();
    // std::cout << "tot M: " << total_magnetization << ", tot_E: " << total_energy << std::endl;
    // std::cout << seed << std::endl;
    // std::cout << std::endl;
    // // best value, forgot seed!!

    return 0;
}

int main()
{   
    int the_magic_seed = 1572032584;
    int n = 20;
    generate_data(the_magic_seed);

    IsingModel q(n);

    return 0;
}