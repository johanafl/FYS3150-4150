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


int compute_and_write_eigendata(double freq, std::string filename, 
                                      bool progress, double rho_min, 
                                      double rho_max, int grid, 
                                      double tol_off_diag, double rho_end, 
                                      double d_rho)
{   
    /* 
    Computing eigenvalues and vectors for matrix approximating the Hamiltonian
    opetrator for a potential well with two electrons interecting with Coulomb
    force.
    NB! Remember to change filename_eig! This program WILL write over existing 
    files.

    Parameters
    ----------
    freq : double
        Oscillator frequency for the harmonic oscillator.

    filename : std::string
        An input for what the filename should be. The final filename will be:
        "eigenvalue_" + filename + ".txt", for the file containing eigenvalues,
        "eigenvector_" + filename + ".txt", for the file containing eigenvectors.

    progress : bool
        If you want to see the progress in the terminal, set this to "true".
        It writes:
        rho_max: xxx rho_end: xxx  n: xxx
        rho_max: xxx rho_end: xxx  n: xxx
        ...
        rho_max: xxx rho_end: xxx  n: xxx

    rho_min : double
        This is where the potential startes. Must be >= 0.

    rho_max : double
        This is where the potential ends for the first iteration.

    grid : int
        This is the number of gridpoints.

    tol_off_diag : double
        The tolleranse for jacobi method.

    rho_end : double
        This is where the iteration stops for infty.

    d_rho : double
        increment for rho_max
    */ 
    std::cout << "Looping over different rho_max values" << std::endl;

    double rho_tmp = rho_max; // for reverting rho_max to original max value
    // Looping over frequencies.
    // Setting up the file.
    std::string filename_eig;
    std::string filename_vec;
    std::ofstream data_file_eig;
    std::ofstream data_file_vec;

    filename_eig = "eigenvalue_" + filename + ".txt";
    filename_vec = "eigenvector_" + filename + ".txt";
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
        double off_diag = -1/(step*step);       // Off-diagonal elements in the matrix we want the eigenvalues of.

        // constructing vector for the diagonal elements            
        arma::vec diag_elements = effective_potential(grid, step, rho_min, freq);

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
            std::cout << "n: " << grid << std::endl;
        }
        
        // Writing ground state eigenvector and value to file.
        data_file_eig << std::setw(6) << rho_max << std::setw(20) << std::setprecision(10) << sorted_diag(0) << std::endl;
        data_file_vec << std::setw(20) << rho_max;
        
        for (int j=0; j<grid; j++)
        {   // writing eigenvector to file
            data_file_vec << std::setw(20) << std::setprecision(10) << ground_state(j);
        }
        
        rho_max = rho_max + d_rho;
        data_file_vec << std::endl;
    }
    rho_max = rho_tmp;
    data_file_eig.close();
    data_file_vec.close();
}


class Eigendata
{
    /*
    All methods contain specific parameters needed for generating good values
    for the different harmonic oscillator frequencies.

    Common values of the class are values constant for all frequencies.

    Runnin any method will append data to an already existing file.
    */
private:
    bool progress = true;        // toggling progress info on/off
    int grid = 200;              // end grid value
    double rho_min = std::pow(10, -7);
    double tol_off_diag = std::pow(10, -5); // tolerance for Jacobi
    double rho_max;
    double rho_end;
    double d_rho;
    double freq;
    std::string filename;

public:
    Eigendata()
    {

    }

    void eigendata_freq_001()
    {
        /*
        Specific parameters for freq = 0.01.

        Unknown good rho max.
        */

        rho_max  = 6.8;     // approximating infinity
        rho_end  = 8;     // end rho_max value for the loop
        d_rho    = 0.2;    // rho step size for incrementing in the loop
        freq     = 0.01; 
        filename = "omega_" + std::to_string(freq);
        compute_and_write_eigendata(freq, filename, progress, rho_min, rho_max,
            grid, tol_off_diag, rho_end, d_rho);
    }

    void eigendata_freq_005()
    {
        /*
        Specific parameters for freq = 0.05.

        Known good rho max.
        */

        rho_max  = 15;     // approximating infinity
        rho_end  = 25;     // end rho_max value for the loop
        d_rho    = 0.2;    // rho step size for incrementing in the loop
        freq     = 0.05; 
        filename = "omega_" + std::to_string(freq);
        compute_and_write_eigendata(freq, filename, progress, rho_min, rho_max,
            grid, tol_off_diag, rho_end, d_rho);
    }

    void eigendata_freq_025()
    {
        /*
        Specific parameters for freq = 0.25.

        Known good rho max.
        */

        rho_max  = 7.7;     // approximating infinity
        rho_end  = 10;     // end rho_max value for the loop
        d_rho    = 0.2;    // rho step size for incrementing in the loop
        freq     = 0.25; 
        filename = "omega_" + std::to_string(freq);
        compute_and_write_eigendata(freq, filename, progress, rho_min, rho_max,
            grid, tol_off_diag, rho_end, d_rho);
    }

    void eigendata_freq_05()
    {
        /*
        Specific parameters for freq = 0.5.

        Unknown good rho max.
        */

        rho_max  = 2;     // approximating infinity
        rho_end  = 8;     // end rho_max value for the loop
        d_rho    = 0.2;    // rho step size for incrementing in the loop
        freq     = 0.5; 
        filename = "omega_" + std::to_string(freq);
        compute_and_write_eigendata(freq, filename, progress, rho_min, rho_max,
            grid, tol_off_diag, rho_end, d_rho);
    }

    void eigendata_freq_1()
    {
        /*
        Specific parameters for freq = 1.

        Unknown good rho max.
        */

        rho_max  = 0.1;     // approximating infinity
        rho_end  = 7;     // end rho_max value for the loop
        d_rho    = 0.2;    // rho step size for incrementing in the loop
        freq     = 1; 
        filename = "omega_" + std::to_string(freq);
        compute_and_write_eigendata(freq, filename, progress, rho_min, rho_max,
            grid, tol_off_diag, rho_end, d_rho);
    }

    void eigendata_freq_5()
    {
        /*
        Specific parameters for freq = 5.

        Unknown good rho max.
        */

        rho_max  = 0.1;     // approximating infinity
        rho_end  = 6;     // end rho_max value for the loop
        d_rho    = 0.2;    // rho step size for incrementing in the loop
        freq     = 5; 
        filename = "omega_" + std::to_string(freq);
        compute_and_write_eigendata(freq, filename, progress, rho_min, rho_max,
            grid, tol_off_diag, rho_end, d_rho);
    }


};




int main()
{    

    Eigendata q;
    q.eigendata_freq_025();
    
    // {0.01, 0.05, 0.25, 0.5, 1, 5}

    return 0;
}