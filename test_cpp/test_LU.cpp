#include <armadillo>
#include <iostream>

int main()
{
    int n_rows = 10;
    int n_cols = 10;
    arma::mat A(n_rows, n_cols);

    A.diag(0)  += 2.0;
    A.diag(1)  += -1.0;
    A.diag(-1) += -1.0;

    arma::mat L;
    arma::mat U;
    arma::lu(L, U, A);

    U.print();
    return 0;
}