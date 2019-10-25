// #include <mpi.h>
#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <armadillo>
double const pi = 3.14159265359; 

double** initial_spin(int n)
{
    time_t seed;
    time(&seed);
    std::mt19937 engine(seed);
    // generating distribution
    std::uniform_int_distribution<int> uniform(0, 1);

    double** spin = new double*[n];

    for (int i=0; i<n; i++)
    {
        spin[i] = new double[n];
        for (int j=0; j<n; j++)
        {
            double spin_val = 2 * uniform(engine) - 1;
            spin[i][j] = spin_val;
        }
    }

    return spin;
}

double total_energy_and_magnetization(double** spin, int n, double& tot_magnet)
{
    /*
    Computing from upper left in matrix and down to the right. 
    */
    double tot_energy = 0;
    
    for (int i=0; i<n-1; i++)
    {
        for (int j=0; j<n-1; j++)
        {
            double spin_now = spin[i][j];
            double spin_right = spin[i][j+1];
            double spin_down = spin[i+1][j];

            tot_energy -= spin_now*(spin_right + spin_down);
            tot_magnet += spin_now;
        }
    }
    
    for (int j=0; j<n-1; j++)
    {
        double spin_now = spin[n-1][j];
        double spin_right = spin[n-1][j+1];
        double spin_down = spin[0][j];

        tot_energy -= spin_now*(spin_right + spin_down);
        tot_magnet += spin_now;
    }

    for (int i=0; i<n-1; i++)
    {
        double spin_now = spin[i][n-1];
        double spin_right = spin[i][0];
        double spin_down = spin[i+1][n-1];

        tot_energy -= spin_now*(spin_right + spin_down);
        tot_magnet += spin_now;
    }

    tot_energy += spin[n-1][n-1]*(spin[n-1][0] + spin[0][n-1]);
    tot_magnet += spin[n-1][n-1];

    return tot_energy;
}

void print_maetrix(double** matrex, int n)
{
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            std::cout << std::setw(4) << matrex[i][j];
        }
        std::cout << std::endl;
    }
}

double enegy_from_flap(double** spin_mat, int row, int col, int n)
{   
    /*This is wrong!!!!!!*/
    double spin_now   = spin[row][col]*(-1);
    // PROV AA DROPPE IF_TEST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!v
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double spin_left  = spin[row][col-1];
    double spin_up    = spin[row-1][col];
    double spin_right = spin[row][col+1];
    double spin_down  = spin[row+1][col];
    if (row == n-1)
    {
        spin_down  = spin[0][col];
    }
    if (row == 0)
    {
        spin_up    = spin[n-1][col];
    }
    if (col == n-1)
    {
        spin_right = spin[row][0];
    }
    if (col == 0)
    {
        spin_left  = spin[row][n-1];
    }

    double delta_energy = spin_now*(spin_up + spin_down + spin_right + spin_left);

    return delta_energy;
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
    double** init_spin = initial_spin(n);


    double matnetization;
    double engery = total_energy_and_magnetization(init_spin, n, matnetization);

    print_maetrix(init_spin, n);

    std::cout << std::endl;
    std::cout << "tot M: " << matnetization << ", tot_E: " << engery << std::endl;

    // HUSK DELETE!!!!!!!!!!
    for (int i=0; i<n; i++)
    {
        delete[] init_spin[i];
    }
    delete[] init_spin;

    return 0;
}
