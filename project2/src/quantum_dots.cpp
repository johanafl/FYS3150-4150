#include "jacobi.h"
#include <fstream>
#include <iomanip>


int main()
{   

    int n = 100;                // grid size
    double rho_min = 0;
    double rho_max = 8.05;      // approx infty
    double step = (rho_max - rho_min)/n;
    double off_diag = -1/(step*step);
    double tol_off_diag = std::pow(10, -5);

    // constructing vector for the diagonal elements
    arma::vec diag(n);
    diag.zeros();

    for (int i = 0; i < n; i++)
    {   // creating the diagonal elements
        diag(i) = 2/(step*step) + (rho_min + i*step)*(rho_min + i*step);
    }

    // constructing tri-diagonal matrix and vector for analytical eigenvalues
    arma::mat A = construct_diag_matrix(n, off_diag, diag);
    arma::vec analytical_eig(n);
    analytical_eig(0) = 3;  // initial value
    
    // using Jacobis method to extract the eigenvalues of the tri-diagonal matrix
    find_eig(n, A, tol_off_diag);
    arma::vec sorted_diag = arma::sort(A.diag(0));


    std::ofstream data_file;
    data_file.open("eigenvalues.txt");
    data_file << "calculated     exact     error\n";

    
    
    for (int i = 0; i < 8; i++)
    {   // writing data to file
        
        analytical_eig(i + 1) = analytical_eig(i) + 4;
        data_file << std::setw(10) << sorted_diag(i);
        data_file << std::setw(10) << analytical_eig(i);
        data_file << std::setw(10) << fabs(analytical_eig(i) - sorted_diag(i));
        data_file << "\n";

    }

    data_file.close();

    return 0;
}