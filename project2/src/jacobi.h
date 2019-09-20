#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>

namespace jac
{
    arma::mat construct_diag_matrix(int n);
    arma::mat construct_diag_matrix(int n, double off_diag, double diag);
    double find_max(int n, arma::mat& A, int& idx_row, int& idx_col);
    void transform(int n, arma::mat& A, arma::mat& R, int idx_row, int idx_col);
    void find_eig(int n, arma::mat& A, double tol_off_diag);
}

