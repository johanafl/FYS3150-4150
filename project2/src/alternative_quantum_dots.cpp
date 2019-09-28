#include "jacobi.h"
// #include <iomanip>
#include <string>

arma::vec generate_diagonals(int grid, double step, double rho_min)
{
    /*
    Creates a vector containing correct diagonal elements for task 2d.

    Parameters
    ----------
    n : int
        Grid size.

    step : double
        Step size.

    rho_min : double
        Minimum rho value.
    */

    arma::vec diag(grid);
    for (int i = 0; i < grid; i++)
    {   // creating the diagonal elements
        diag(i) =   2/(step*step) + (rho_min + i*step)*(rho_min + i*step);
    }
    return diag;
}

int main()
{   
    /* 
    Calculation variables. 
    NB! Must be specified correctly before program running the program! 
    NB! Remember to change filename! This program WILL write over existing 
    files.
    */ 
    bool progress        = true;  // boolean for toggling progress info on/off
    int grid             = 100;   // grid size
    double rho_min       = std::pow(10, -7);
    double rho_max       = 5;     // approximating infinity
    double tol_off_diag  = std::pow(10, -10); // tolerance for Jacobi
    std::string filename = "eigenvalues2.txt";

    // loop-specific values
    double rho_tmp = rho_max; // for reverting rho_max to original max value
    double rho_end = 5.6;     // end rho_max value for the loop
    int grid_tmp   = grid;    // for reverting n to original max value
    int grid_end   = 130;     // end grid value
    double d_rho   = 0.1;     // rho step size for incrementing in the loop
    int d_grid     = 10;      // grid step size
    int num_eig    = 8;       // number of eigenvalues to write to file per rho

    // Setting up the file
    std::ofstream data_file;
    data_file.open(filename, std::ios_base::app);
    data_file << std::setw(10) << "n | rho ->";

    while (rho_max <= rho_end)
    {
        data_file << std::setw(20) << rho_max;
        rho_max += d_rho;
    }
    data_file << std::endl;
    rho_max = rho_tmp;

    // Doing the calculations
    while (grid <= grid_end)
    {
        // Looping over grid values.
        data_file << std::setw(10) << grid;

        while (rho_max <= rho_end)
        {
            // Looping over maximum values of rho (approximating infinity).
            double eig = 3;    // initial analytical eigenvalue
            double step = (rho_max - rho_min)/grid; // Step size
            double off_diag = -1/(step*step);       // Off-diagonal elements in the matrix we want the eigenvalues of.

            // constructing vector for the diagonal elements            
            arma::vec diag = generate_diagonals(grid, step, rho_min);

            // constructing tri-diagonal matrix
            arma::mat A = construct_diag_matrix(grid, off_diag, diag);
            
            // using Jacobis method to extract the eigenvalues of the tri-diagonal matrix
            arma::mat R = find_eig(grid, A, tol_off_diag);
            // The eigenvalues are not nessesarely sorted, and we want the ground state.
            arma::vec sorted_diag = arma::sort(A.diag(0));

            if (progress)
            {   // for printing progress data
                std::cout << "rho_max: " << rho_max << " rho_end: ";
                std::cout << rho_end << "  ";
                std::cout << "n: " << grid << " n_end: ";
                std::cout << grid_end << std::endl;
            }

            arma::vec error_rho_n(num_eig);
            for (int i = 0; i < num_eig; i++)
            {   
                // writing data to file
                error_rho_n(i) = fabs(eig - sorted_diag(i));
                eig = eig + 4;
            }
            data_file << std::setw(20) << std::setprecision(10) << error_rho_n.max();
            rho_max = rho_max + d_rho;
        }
        grid += d_grid;
        data_file << std::endl;
    }
    data_file.close();

    return 0;
}