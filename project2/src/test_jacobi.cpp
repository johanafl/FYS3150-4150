#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "jacobi.h"


TEST_CASE("test_find_max")
{
    /*
    Test function that tests the most important functionality of the Jacobi class.

    This test works.
    */
    
    
    int n = 5;                  // dimension of test matrix
    double maximum = 5555;      // arbitrary max value in matrix
    
    arma::mat A(n, n);
    A.zeros();
    
    for (int i = 0; i < 5; i++)
    {   // loops over rows
        
        for (int j = 0; j < 5; j++)
        {   // loops over columns
            A(i, j) = i*j*1.5;  // populating the array with arbitrary values
        }
    }

    A(3, 2) = maximum;  // placing the max value
    
    int idx_row;
    int idx_column;
    double max_val = find_max(n, A, idx_row, idx_column);

    REQUIRE(max_val == maximum);
    REQUIRE(idx_row == 3);
    REQUIRE(idx_column == 2);

}


TEST_CASE("test_inner_product_conserved")
{
    /*
    Checks that the inner product is preserved after rotation with the
    transformation matrix. Sets up an anti diagonal matrix where the columns are
    R^(4x4) basis vectors. Checks that the resulting inner product of all
    combinations of different vectors is 0.

    This test works.
    */
    
    int n = 4;          // dimension of matrix
    int idx_row = 0;
    int idx_col = 3;
    double tol = std::pow(10, -10);
    
    arma::mat R(n, n);      // matrix for storing the eigenvectors of A
    R.zeros();
    R.diag(0) += 1;
    
    arma::mat A(n, n);
    A.zeros();
    A(0, 3) = 1;
    A(1, 2) = 1;
    A(2, 1) = 1;
    A(3, 0) = 1;
    

    transform(n, A, R, idx_row, idx_col);

    for (int i = 0; i < n; i++)
    {   // loops over all column vectors
        
        for (int j = 0; j < n; j++)
        {   // loops over all column vectors
            
            if (i != j)
            {   // discards equal vectors
                REQUIRE(fabs(arma::dot(A.col(i), A.col(j))) < tol);
            }
        }
    }
}


TEST_CASE("test_find_eig")
{
    /*
    Checks that the find_eig function finds eigenvalues up to a fixed tolerance.
    */
    
    double pi = 3.14159265358979323846;
    int n = 5;                                  // dimension of matrix
    double step     = 1.0/n;
    double diag     = 2/(step*step);            // diagonal elements
    double off_diag = -1/(step*step);           // off-diagonal elements
    double tol_off_diag = std::pow(10, -15);    // tolerance of off-diagonal elements
    double tol_eig  = std::pow(10, -4);         // tolerance of eigenvalues
    
    arma::mat A = construct_diag_matrix(n, off_diag, diag);
    arma::vec eigenvalues(n);

    
    for (int i = 1; i <= n; i++)
    {   // generating analytical eigenvalues
        eigenvalues(i - 1) = diag + 2*off_diag*std::cos(i*pi/(n + 1));
    }
    
    arma::mat R = find_eig(n, A, tol_off_diag); // numerical eigenvalues

    arma::vec sorted_diag = arma::sort(A.diag(0));
    
    for (int i = 0; i < n; i++)
    {   // checking that analytical and numerical results match to a given tolerance
        
        REQUIRE(fabs(eigenvalues(i) - sorted_diag(i)) < tol_eig);
    }
}
