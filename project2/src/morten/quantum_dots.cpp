// #include <iostream>
#include "jacobi.h"



int main()
{   

    int n = 10;
    double rho_min = 0;
    double rho_max = 5;    // 5 approx infty
    double step = (rho_max - rho_min)/n;
    double off_diag = -1/(step*step);

    arma::vec diag(n);
    diag.zeros();
    diag.print();

        
    for (int i = 0; i < n; i++)
    {   // creating the diagonal elements
        diag(i) = 2/(step*step) + (rho_min + i*step)*(rho_min + i*step);
    }

    arma::mat A = jac::construct_diag_matrix(n, off_diag, diag);
    jac::hei();

    std::cout << "lol" << std::endl;
    return 0;
}
