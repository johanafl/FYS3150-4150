#ifndef JACOBI_H
#define JACOBI_H

#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>


arma::mat construct_diag_matrix(int n);
arma::mat construct_diag_matrix(int n, double off_diag, double diag);
arma::mat construct_diag_matrix(int n, double off_diag, arma::vec diag);
double find_max(int n, arma::mat& A, int& idx_row, int& idx_col);
void transform(int n, arma::mat& A, arma::mat& R, int idx_row, int idx_col);
arma::mat find_eig(int n, arma::mat& A, double tol_off_diag);
void test_function_for_checking_header_implementation();

const double pi = 3.14159265358979323846;

#endif