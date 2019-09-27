#include "jacobi.h"
#include <fstream>
#include <iomanip>


void quantum_data()
{   /*

    Loops over rho_max and n (grid) values, calculates eigenvalues by the Jacobi
    method, and writes data to a text file for visualization with Python.
    */

    bool progress = true;   // boolean for toggling progress info on/off
    double eig;             // analytical eigenvalue
    int n = 100;            // grid size
    double rho_min = 0;
    double rho_max = 5;         // approx infty
    double tmp = rho_max;
    double tol_off_diag = std::pow(10, -5); // tolerance for Jacobi
    double step;
    double off_diag;
    

    // loop-specific values
    double d_rho   = 0.1;         // rho step size for incrementing in the loop
    int num_eig    = 8;           // number of eigenvalues to write to file
    double rho_end = 10;        // end rho_max value for the loop
    double num_rho = (rho_end - rho_max)/d_rho; // number of rho_max values tested
    int n_end = 410;    // end grid value
    int dn    = 10;     // grid step size
    int num_n = (n_end - n)/dn;


    // generating data file
    std::ofstream data_file;
    data_file.open("eigenvalues.txt", std::ios_base::app);
    data_file << num_eig << " " << num_rho << " " << num_n << "\n";
    data_file << std::setw(20) << "calculated";
    data_file << std::setw(20) << "exact";
    data_file << std::setw(20) << "error";
    data_file << std::setw(20) << "rho_max";
    data_file << std::setw(20) << "n\n";

    if (progress)
    {   // for printing progress data
        std::cout << "calculating" << std::endl;
        std::cout << "===========" << std::endl;
    }

    while (n < n_end)
    {
        while (rho_max < rho_end)
        {   // looping over rho_max values and writing data to file
            
            eig = 3;    // initial analytical eigenvalue
            step = (rho_max - rho_min)/n;
            off_diag = -1/(step*step);

            // constructing vector for the diagonal elements
            arma::vec diag(n);
            diag.zeros();

            for (int i = 0; i < n; i++)
            {   // creating the diagonal elements
                diag(i) = 2/(step*step) + (rho_min + i*step)*(rho_min + i*step);
            }

            // constructing tri-diagonal matrix
            arma::mat A = construct_diag_matrix(n, off_diag, diag);
            
            // using Jacobis method to extract the eigenvalues of the tri-diagonal matrix
            find_eig(n, A, tol_off_diag);
            arma::vec sorted_diag = arma::sort(A.diag(0));

            if (progress)
            {   // for printing progress data
                std::cout << "rho_max: " << rho_max << " rho_end: ";
                std::cout << rho_end << "  ";
                std::cout << "n: " << n << " n_end: ";
                std::cout << n_end << std::endl;
            }


            for (int i = 0; i < num_eig; i++)
            {   // writing data to file
                
                data_file << std::setw(20) << std::setprecision(10) << sorted_diag(i);
                data_file << std::setw(20) << std::setprecision(10) << eig;
                data_file << std::setw(20) << std::setprecision(10) << fabs(eig - sorted_diag(i));
                data_file << std::setw(20) << rho_max;
                data_file << std::setw(20) << n;
                data_file << "\n";
                eig = eig + 4;

            }


            rho_max = rho_max + d_rho;

        }
        
        rho_max = 5;
        n += dn;
    }
    
    data_file.close();

}

int main()
{   
    quantum_data();
    return 0;
}