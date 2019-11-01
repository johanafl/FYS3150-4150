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
    double spin_here = spin(row, col, true);
    double delta_energy = 2*spin_here*(spin(row-1, col, true) + spin(row+1, col, true)
                            + spin(row, col+1, true) + spin(row, col-1, true));

    
    if (metropolis_random < exp_delta_energy[(int) (delta_energy + 8)])
    {   // checks if energy difference is positive and the metropolis condition true
        spin(row, col, true) *= -1;
        total_energy         += delta_energy;
        total_magnetization  += -2*spin_here;
    }
}


class IsingModel
{

private:
    int n;              // matrix is of dimension nxn
    int mc_iterations;  // number of Monte Carlo iterations
    double* exp_delta_energy = new double[17];  // pre-calculated energy values
    double J = 1;

    int row;    // row index, will be randomly drawn
    int col;    // column index, will be randomly drawn
    double metropolis_random;   // metropolis condition, will be randomly drawn
    double total_energy;
    double total_magnetization;

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

    // defining the spin matrix
    CircularMatrix spin;

    void mc_iteration_convergence(double temp)
    {   /*
        Runs the spin flip a given amount of times. Generates data for finding
        how many iterations is needed for convergence. Keeps the energy values
        without averaging. Only runs for a single temperature value.
        */

        for (int j = 0; j < mc_iterations; j++)
        {   // loops over n*n spin flips a given amount of times
            // saves relevant data for each iteration

            iterate_spin_flip(temp);
            E_convergence_data << std::setw(15) << total_energy;
            M_convergence_data << std::setw(15) << total_magnetization;
        }
    }

    void mc_iteration_stable(double temp)
    {   /*
        Runs the spin flip a given amount of times. Does not keep every energy
        value, but calculates the average.

        Parameters
        ----------
        temp : double
            Temperature value.
        */

        sum_total_energy = 0;
        sum_total_energy_squared = 0;
        sum_total_magnetization  = 0;
        sum_total_magnetization_absolute = 0;
        sum_total_magnetization_squared  = 0;

        for (int j = 0; j < mc_iterations; j++)
        {   // loops over n*n spin flips a given amount of times

            iterate_spin_flip(temp);

            sum_total_energy += total_energy;
            sum_total_energy_squared += total_energy*total_energy;
            sum_total_magnetization  += total_magnetization;
            sum_total_magnetization_absolute += std::fabs(total_magnetization);
            sum_total_magnetization_squared  += total_magnetization*total_magnetization;

        }
        
        sum_total_energy /= mc_iterations;
        sum_total_energy_squared /= mc_iterations;
        sum_total_magnetization  /= mc_iterations;
        sum_total_magnetization_absolute /= mc_iterations;
        sum_total_magnetization_squared  /= mc_iterations;

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

public:    
    IsingModel(int spin_mat_dim, int mc_iterations_input, long seed) : uniform_discrete(0, spin_mat_dim - 1), engine(seed), uniform_continuous(0, 1), spin(spin_mat_dim, seed)
    {   /*
        Parameters
        ----------
        spin_mat_dim : int
            The spin matrix is of dimenstion spin_mat_dim x spin_mat_dim.

        mc_iterations_input : int
            The number of Monte Carlo iterations.

        seed : long
            Seed for the PRNG.
        */
        
        n = spin_mat_dim;
        mc_iterations = mc_iterations_input;
        total_energy_and_magnetization(spin, n, total_energy, total_magnetization);
        
        // initialising data files
        E_convergence_data.open("data_files/E_convergence_data.txt", std::ios_base::app);
        M_convergence_data.open("data_files/M_convergence_data.txt", std::ios_base::app);
        ising_model_data.open("data_files/ising_model_data.txt", std::ios_base::app);
    }

    void iterate_temperature(double initial_temp, double final_temp,
        double dtemp, bool convergence)
    {   /*
        Iterates over a given set of temperature values.

        Parameters
        ----------
        initial_temp : double
            Start temperature value.

        final_temp : double
            End temperature value.

        dtemp : double
            Temperature step length.

        convergence : bool
            This class generates data, either for the convergence phase or for
            the stable phase. convergence toggles this.
        */

        if (convergence)
        {
            E_convergence_data << "mc_iterations: " << mc_iterations;
            E_convergence_data << " spin_matrix_dim: " << n;
            E_convergence_data << std::endl;
            M_convergence_data << "mc_iterations: " << mc_iterations;
            M_convergence_data << " spin_matrix_dim: " << n;
            M_convergence_data << std::endl;
        }
    
        for (double temp = initial_temp; temp <= final_temp; temp += dtemp)
        {   // looping over temperature values

            // pre-calculated exponential values
            exp_delta_energy[0]  = std::exp(8*J/temp);
            exp_delta_energy[4]  = std::exp(4*J/temp);
            exp_delta_energy[8]  = 1;
            exp_delta_energy[12] = std::exp(-4*J/temp);
            exp_delta_energy[16] = std::exp(-8*J/temp);

            if (convergence)
            {   // generates and writes convergence data
                
                // writing temperature values in the first column
                E_convergence_data << std::setw(15) << temp;
                M_convergence_data << std::setw(15) << temp;
                
                mc_iteration_convergence(temp);
                
                E_convergence_data << std::endl;
                M_convergence_data << std::endl;
            }
            else
            {   // generates and writes stable data
                
                mc_iteration_stable(temp);
                ising_model_data << std::setw(20) << std::setprecision(15) << temp;
                ising_model_data << std::setw(20) << std::setprecision(15) << sum_total_energy;
                ising_model_data << std::setw(20) << std::setprecision(15) << sum_total_energy_squared;
                ising_model_data << std::setw(20) << std::setprecision(15) << sum_total_magnetization;
                ising_model_data << std::setw(20) << std::setprecision(15) << sum_total_magnetization_squared;
                ising_model_data << std::setw(20) << std::setprecision(15) << sum_total_magnetization_absolute;
                ising_model_data << std::endl;
            }     
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




    return 0;
}

int main()
{   
    int the_magic_seed = 1572032584;
    int spin_matrix_dim = 2;
    int mc_iterations = 1e3;
    bool convergence = true;
    
    double initial_temp = 1;
    double final_temp = 1;
    double dtemp = 1;       // temperature step length

    time_t seed;
    time(&seed);
    
    IsingModel q(spin_matrix_dim, mc_iterations, seed);
    q.iterate_temperature(initial_temp, final_temp, dtemp, convergence);

    return 0;
}