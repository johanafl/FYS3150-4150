// #include <mpi.h>
#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <armadillo>
double const pi = 3.14159265359;



double** initial_spin(int n)
{   /*
    Creates a n x n matrix where last row and last column are copies of the
    first row and first column (padding).

    Parameters
    ----------
    n : int
        Grid dimension.

    Returns
    -------
    spin : double**
        Matrix of random spin values (+- 1).
    */
    
    int seed = 1337;
    // time_t seed;
    // time(&seed);
    std::mt19937 engine(seed);
    std::uniform_int_distribution<int> uniform(0, 1);

    double** spin = new double*[n+2];

    for (int i = 0; i < n + 2; i++)
    {   // creates an array for every element in the spin array
        // n + 1 here to create the last column, even though we overwrite the
        // values later
        
        spin[i] = new double[n+2];
        
        for (int j = 1; j < n + 1; j++)
        {   // creates +- 1 values and inserts into array
            
            double spin_val = 2*uniform(engine) - 1;  // +- 1
            spin[i][j] = spin_val;
        }
    }

    std::cout << &spin[1][1] << std::endl;

    // for (int i = 0; i < n; i++)
    // {   // copies first row and first column into the last row and last
    //     // i < n instead of i < n + 1 drops the last unused element.
        
    //     &spin[n][i] = &spin[0][i];
    //     &spin[i][n] = &spin[i][0];
    // }

    return spin;
}

double total_energy_and_magnetization(double** spin, int n, double& tot_magnet)
{   /*
    Computing from upper left in matrix and down to the right.

    Parameters
    ----------
    spin : double**
        n x n matrix represented by a double pointer.
    
    n : int
        Grid dimension.

    tot_magnet : double refernce
        Total magnetic moment.
    */

    double tot_energy = 0;
    
    for (int i = 0; i < n; i++)
    {   // looping over rows
        
        for (int j = 0; j < n; j++)
        {   // looping over columns

            double spin_here  = spin[i][j];     // spin of the current indices
            double spin_below = spin[i+1][j];   // below
            double spin_right = spin[i][j+1];   // right

            tot_energy -= spin_here*(spin_right + spin_below);
            tot_magnet += spin_here;
        }
    }
    
 
    return tot_energy;
}

void print_matrix(double** matrix, int n, bool pad)
{   /*
    Prints matrix.

    Parameters
    ----------
    matrix : double**
        n x n matrix represented by a double pointer (n+1 x n+1 with pad).

    n : int
        Grid dimension.

    pad : bool
        Option to toggle print of padding on/off.
    */

    int m = n;

    if (pad)
    {   // for the extra row and column
        m = n + 1;
    }

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            std::cout << std::setw(4) << matrix[i][j];
        }
        std::cout << std::endl;
    }
}

void metropolis_flap(double** spin, double& tot_energy, double& tot_magnet,
        int row, int col, double metro_con, double T)
{   /*
    
    Parameters
    ----------
    spin : double**
        Matrix of spin values.
    
    tot_energy : double reference
        Total energy of all the spins.

    tot_magnet : double reference
        Total magnetic moment of the system.

    row : int
        A randomly chosen row index.

    col : int
        A randomly chosen column index.

    metro_con : double
        Random variable drawn from a uniform distribution on the interval
        [0, 1). Metropolis condition.
    */

    std::cout << "row: " << row << " col: " << col << std::endl;
    
    double spin_here  = -1*spin[row][col];  // flippity flap
    double spin_left  = spin[row][col-1];
    double spin_above = spin[row-1][col];
    double spin_right = spin[row][col+1];
    double spin_below = spin[row+1][col];

    double delta_energy = 2*spin_here*(spin_above + spin_below + spin_right + spin_left);

    // if (delta_energy <= 0)
    // {   // accept new energy if difference is negative
        
    //     spin[row][col] *= -1;
    //     tot_energy += delta_energy;
    //     tot_magnet += spin_here;
    // }
    
    // else if ( (delta_energy > 0) and (metro_con < std::exp(delta_energy/T)) )
    // {   // checks if energy difference is positive and the metropolis condition true
        
    //     spin[row][col] *= -1;
    //     tot_energy += delta_energy;
    //     tot_magnet += spin_here;
    // }
    

}


// void mc_integration(arma::mat spin_mat)
// {
//     double sum_energy = 0;
//     for (int i=0; i<n; i++)
//     {

//     }
// }



int main()
{   
    // PARALLELL!!!!!!
    //
    // MPI_Init(NULL, NULL);
    // int world_rank;
    // int world_size;
    // MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // MPI_Finalize();

    int n = 3;
    double temperature = 1;
    
    int seed = 1337;
    // time_t seed;
    // time(&seed);
    std::mt19937 engine(seed);
    std::uniform_int_distribution<int> uniform_discrete(0, n - 1);
    std::uniform_real_distribution<double> uniform_continuous(0, 1);




    double** init_spin = initial_spin(n);

    // double magnetization;
    // double energy = total_energy_and_magnetization(init_spin, n, magnetization);

    // print_matrix(init_spin, n, true);

    // std::cout << std::endl;
    // std::cout << "tot M: " << magnetization << ", tot_E: " << energy << "\n" << std::endl;
    
    // for (int i = 0; i < 10; i++)
    // {   
    //     int row = uniform_discrete(engine);
    //     int col = uniform_discrete(engine);
    //     double metro_cond = uniform_continuous(engine);
        
    //     metropolis_flap(init_spin, energy, magnetization, row, col, metro_cond, temperature);
    // }


    // print_matrix(init_spin, n, true);

    // std::cout << std::endl;
    // std::cout << "tot M: " << magnetization << ", tot_E: " << energy << std::endl;

    // HUSK DELETE!!!!!!!!!!
    for (int i = 0; i < n + 1; i++)
    {
        delete[] init_spin[i];
    }
    delete[] init_spin;


    return 0;
}
