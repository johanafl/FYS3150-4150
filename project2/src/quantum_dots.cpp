#include "jacobi.h"


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
    
    // using Jacobis method to extract the eigenvalues of the tri-diagonal matrix
    find_eig(n, A, tol_off_diag);
    arma::vec sorted_diag = arma::sort(A.diag(0));

    for (int i = 0; i < 5; i++)
    {   
        analytical_eig(i) = diag(i) + 2*off_diag*std::cos(i*3.1415/(n + 1));
        
        std::cout << sorted_diag(i) << "  " << analytical_eig(i) << std::endl;
    }

    return 0;
}