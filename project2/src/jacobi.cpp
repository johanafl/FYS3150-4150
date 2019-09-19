#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <chrono>


arma::mat construct_diag_matrix(int n) {
    /*
    Function for creating a tridiagonal matrix. Second derivative discretized.
    */

    arma::mat A(n, n);
    A.zeros();

    A.diag(0)  += 2;
    A.diag(-1) += -1.0;
    A.diag(1)  += -1.0;

    return A;
}

arma::mat construct_diag_matrix(int n, double up_low_diag, arma::vec diag) {
    /*
    Function for creating a diagonal armadillo matrix.
    */
    arma::mat A(n, n);
    A.zeros();

    A.diag(0) = diag;
    A.diag(-1) += up_low_diag;
    A.diag(1) += up_low_diag;

    return A;
}

double find_max(int n, arma::mat A, int* idx_row, int* idx_column)
{   
    /*
    Finds the max value of the off-diagonal elements of the matrix. Saves
    the indices of the max value.

    Parameters
    ----------
    A : arma::mat
        Matrix of which we want to find the maximum element.

    idx_row : int
        Index of the row with maximum element.

    idx_column : int
        Index of the column with maximum element.


    Returns
    -------
    max_val
    */

    double max_val = 0;     // final max value
    double matrix_elemet;   // placeholder for matrix element value to check

    for (int row=0; row<n; row++)
    {   // loops over rows
        
        for (int column=0; column<n; column++)
        {   // loops over columns
            
            if (row != column)
            {   // exclude diagonal elements
                matrix_elemet = fabs(A(row,column));
                
                if (matrix_elemet > max_val)
                {
                    max_val = matrix_elemet;
                    * idx_row = row;
                    * idx_column = column;
                }
            }
        }
    }
    
    return max_val;
}

void transform(int n, arma::mat& A, int idx_row, int idx_col) {
    /*
    This function rotatetes the (symmetric) matrix A by an angle theta in the 
    plane spanned by the unit vectors e_l and e_k, where 
    e_i = [0, 0, ..., 0, 1, 0, ..., 0]. The angle theta is defined s.t. the 
    elements a_lk = a_kl = 0 after the transformation.

    Parameters
    ----------
    A : arma::mat&
        Reference to matrix A

    idx_row : int
        The index k of one of the unit vectors that are used in the rotation.

    idx_col : int
        The index l of the other unit vector that are used in the rotation.

    */


    double a_kk = A(idx_row, idx_row);
    double a_ll = A(idx_col, idx_col);
    double a_lk = A(idx_col, idx_row);

    double tau = (a_ll - a_kk)/(2*a_lk); // We have defined tau = cot(2*theta), where theta is unknown. We are choosing thata s.t. the element b_kl = b_lk = 0, where B is the new matrix after the transformation with elements b_ij.
    double t;
    
    // to counteract numerical precision errors when tau is large
    if (tau >= 0)
    {
        t = 1.0/(tau + std::sqrt(1 + tau*tau));
    }
    else
    {
        t = -1.0/(-tau + std::sqrt(1 + tau*tau));
    }
    
    double c   = 1/std::sqrt(1 + t*t);
    double s   = t*c;
    double a_ki;

    for (int i=0; i<n; i++)
    {   // elements kk, ll, lk, kl should not be changed, but it is fixed after
        // the loop instead of using an if-test in the loop (vectorization)
        
        a_ki = A(idx_row, i);    // since this matrix element will be overwritten
        
        A(i, idx_row) = A(idx_row, i) = a_ki*c - A(idx_col, i)*s;
        A(i, idx_col) = A(idx_col, i) = A(idx_col, i)*c + a_ki*s;
    }

    A(idx_col, idx_row) = A(idx_row, idx_col) = 0; // by virtue of the algorithm (this is how we found theta!).
    A(idx_row, idx_row) = a_kk*c*c - 2*a_lk*c*s + a_ll*s*s;
    A(idx_col, idx_col) = a_ll*c*c + 2*a_lk*c*s + a_kk*s*s; // check with Johan if + is correct

}


arma::mat find_eig(int n, arma::mat A, double tol)
{
    /*
    Not implemented yet.
    */
    
    arma::mat B = A;
    
    int idx_col;
    int idx_row;
    double max_val = find_max(n, B, & idx_col, & idx_row);
    
    while (max_val > tol)
    {
        transform(n, B, idx_col, idx_row);    // no B is returned, runs forever
        max_val = find_max(n, B, & idx_col, & idx_row);
    }
    return B;
}

void test_find_max() {
    /*
    Test function that tests the most important functionality of the Jacobi class.

    This test works.
    */
    
    
    int n = 5;                  // dimension of test matrix
    double maximum = 5555;      // arbitrary max value in matrix
    
    arma::mat A(n, n);
    A.zeros();
    
    for (int i=0; i<5; i++)
    {   // loops over rows
        for (int j=0; j<5; j++)
        {   // loops over columns
            A(i, j) = i*j*1.5;  // populating the array with arbitrary values
        }
    }

    A(3, 2) = maximum;  // placing the max value
    
    int idx_row;
    int idx_column;
    double max_val = find_max(n, A, & idx_row, & idx_column);
    
    if (max_val != maximum)
    {
        A.print();
        std::cout << "Wrong element! Maximum value was " << maximum << ". Got "
                  << max_val << std::endl;
    }
    if ((idx_row != 3) || (idx_column != 2))
    {
        A.print();
        std::cout << "Wrong indecies! Indecies was " << 3 << " and " << 2 
                  << ". Got " << idx_row << " and " << idx_column << std::endl;
    }


    std::cout << max_val << std::endl;  // debug, remove later
}


void test_inner_product_conserved() {
    /*
    Checks that the inner product is preserved after one rotation with the
    transformation matrix.
    */
    
    int n = 3;  // dimension of matrix
    
    arma::mat A(n, n);
    A.zeros();
    
    A(0, 2) = 1;
    A(1, 1) = 1;
    A(2, 0) = 1;
    
    double tol  = 1;
    int idx_row = 0;
    int idx_col = 2;

    transform(n, A, idx_row, idx_col);

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            std::cout << arma::dot(A.col(i), A.col(j)) << std::endl;
        }
    }

    // std::cout << inner_prod << std::endl;
    
    // if (fabs(inner_prod) > tol)
    // {
    //     std::cout << "Something went wrong!" << std::endl;
    // }
}

void test_find_eig() {
    /*
    Checks that the find_eig function finds eigenvalues up to a set error eps.
    */
    
    double tol_eig  = 0.01; // tolerance of non-diagonal elements
    double tol_test = 1;    // tolerance of eigenvalues
    int n = 5;              // dimension of matrix
    
    arma::mat A = construct_diag_matrix(n);
    arma::vec eigenvalues(n);
    
    double d = 2;   // diagonal elements
    double a = -1;  // off-diagonal elements
    double pi = 3.14159265358979323846;
    
    for (int i=0; i<n; i++)
    {   // generating analytical eigenvalues
        eigenvalues(i) = d + 2*a*std::cos(i*pi/(n+1));
    }
    
    arma::mat B = find_eig(n, A, tol_eig);
    
    for (int i=0; i<n; i++)
    {   // checking that analytical and numerical results match to a given tolerance
        if (fabs(eigenvalues(i) - B(i,i)) > tol_test)
        {
            std::cout << "Something went wrong!\n" << "Exact eigenvalue = " 
                      << eigenvalues(i) << ", computed eigenvalue = " << B(i,i)
                      << "." << std::endl;
        }
    }
}


int main(int argc, char* argv[]) {
    // int n = 5;
    // double stepsize = 1/(n+1);


    // arma::mat A = construct_diag_matrix(n);
    // A.print();
    // test_find_max();
    test_inner_product_conserved();
    // test_find_eig();


    return 1;
}