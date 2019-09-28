#include "jacobi.h"
#include <iomanip>
#include <string>

arma::vec effective_potential(int grid, double step, double rho_min, double freq)
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

    freq : double
        A freqency in front of the harmonic oscillator potential. 
        For the case of two electrons this normalizes the constants in the 
        coulomb potential and makes it unitless.
    */

    arma::vec diag_elements(grid);
    for (int i = 0; i < grid; i++)
    {   
        // creating the diagonal elements
        diag_elements(i) =   2/(step*step) 
                  + freq*freq*(rho_min + i*step)*(rho_min + i*step) // Harmonic oscillator potential (with frequency)
                  + 1/(rho_min + i*step);                     // Coulomb potential
    }

    return diag_elements;
}


int main()
{   
    /* 
    Calculation variables. 
    NB! Must be specified correctly before program running the program! 
    NB! Remember to change filename_eig! This program WILL write over existing 
    files.
    */ 
    bool progress       = true;  // boolean for toggling progress info on/off
    double rho_min      = std::pow(10, -7);
    double rho_max      = 5;     // approximating infinity
    int grid        = 100;   // end grid value
    double tol_off_diag = std::pow(10, -5); // tolerance for Jacobi

    // loop-specific values
    double rho_tmp = rho_max; // for reverting rho_max to original max value
    double rho_end = 5.2;     // end rho_max value for the loop
    double d_rho   = 0.1;     // rho step size for incrementing in the loop


    // Doing the calculations
    std::string filename_eig;
    std::string filename_vec;
    std::ofstream data_file_eig;
    std::ofstream data_file_vec;
    // These are the harmonic oscillator frequencies we will look at. For 
    // omega = 0.25 and omega = 0.05, we have analytical solutions for 
    // eigenvector and value of the ground state.
    double freq[6] = {0.01, 0.05, 0.25, 0.5, 1, 5}; 

    // progress information
    std::cout << "looping over frequency and rho max" << std::endl;

    for (int i = 0; i < 6; i++)
    {   
        // Looping over frequencies.
        // Setting up the file.
        filename_eig = "eigenvalue_omega_" + std::to_string(freq[i]) + ".txt";
        filename_vec = "eigenvector_omega_" + std::to_string(freq[i]) + ".txt";
        data_file_eig.open(filename_eig, std::ios_base::app);
        data_file_vec.open(filename_vec, std::ios_base::app);

        data_file_eig << std::setw(6) << "rho" << std::setw(20) << "eigenvalue" << std::endl;
        data_file_vec << std::setw(20) << "rho | eigenvector ->" << std::endl;

        while (rho_max <= rho_end)
        {
            /* 
            Looping over interesting maximum values for rho (approximation of infinity).
            */
            double step = (rho_max - rho_min)/grid; // Step size
            double off_diag = -1/(step*step);           // Off-diagonal elements in the matrix we want the eigenvalues of.

            // constructing vector for the diagonal elements            
            arma::vec diag_elements = effective_potential(grid, step, rho_min, freq[i]);

            // constructing tri-diagonal matrix
            arma::mat diag_matrix = construct_diag_matrix(grid, off_diag, diag_elements);
            
            // using Jacobis method to extract the eigenvalues of the tri-diagonal matrix
            arma::mat eigenvectors = find_eig(grid, diag_matrix, tol_off_diag);

            // The eigenvalues are not necessarily sorted, and we want the ground state.
            arma::uvec sort_indices = arma::sort_index(diag_matrix.diag(0));
            arma::vec sorted_diag = arma::sort(diag_matrix.diag(0));
            arma::vec ground_state = eigenvectors.col(sort_indices(0));

            if (progress)
            {   // for printing progress data
                std::cout << "rho_max: " << rho_max << " rho_end: ";
                std::cout << rho_end << "  ";
                std::cout << "n: " << grid;
                std::cout << ", omega: " << freq[i] << std::endl;
            }
            
            // Writing ground state eigenvector and value to file.
            data_file_eig << std::setw(6) << rho_max << std::setw(20) << std::setprecision(10) << sorted_diag(0) << std::endl;
            data_file_vec << std::setw(20) << rho_max;
            for (int j=0; j<grid; j++)
            {
                data_file_vec << std::setw(20) << std::setprecision(10) << ground_state(j);
            }
            rho_max = rho_max + d_rho;
            data_file_vec << std::endl;
        }
        rho_max = rho_tmp;
        data_file_eig.close();
        data_file_vec.close();
    }
    
    return 0;
}