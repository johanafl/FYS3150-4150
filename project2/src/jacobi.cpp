#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <chrono>


arma::mat construct_diag_matrix(int n)
{
    /*
    Function for creating a tridiagonal matrix. Second derivative discretized.

    Parameters
    ----------
    n : int
        Dimension of matrix.

    Returns
    -------
    A : arma::mat
        Tridiagonal matrix representing the second derivative discretized.
    */

    arma::mat A(n, n);
    A.zeros();

    A.diag(0)  += 2;
    A.diag(-1) += -1.0;
    A.diag(1)  += -1.0;

    return A;
}


arma::mat construct_diag_matrix(int n, double off_diag, double diag)
{
    /*
    Function for creating a tridiagonal matrix. With vector inputs for
    populating the upper and lower diagonal.

    Parameters
    ----------
    n : int
        Dimension of matrix.

    off_diag : double
        Value for the upper and lower diagonal elements.

    diag : double
        Value for the diagonal elements.
    */

    arma::mat A(n, n);
    
    A.zeros();
    A.diag(0)  += diag;
    A.diag(-1) += off_diag;
    A.diag(1)  += off_diag;

    return A;
}


double find_max(int n, arma::mat& A, int& idx_row, int& idx_col)
{   
    /*
    Finds the max value of the off-diagonal elements of the matrix. Saves
    the indices of the max value.

    Parameters
    ----------
    n : int
        Dimension of matrix.

    A : arma::mat&
        Reference to matrix of which we want to find the maximum value element.

    idx_row : int&
        Reference to row index value where the row index of the maximum value
        will be stored.

    idx_col : int&
        Reference to column index value where the column index of the maximum
        value will be stored.


    Returns
    -------
    max_val : double
        The maximum value in matrix A.
    */

    double max_val = 0;     // final max value
    double matrix_elemet;   // placeholder for matrix element value to check

    for (int row = 0; row < n; row++)
    {   // loops over rows
        
        for (int col = 0; col < n; col++)
        {   // loops over columns
            
            if (row != col)
            {   // exclude diagonal elements
                matrix_elemet = fabs(A(row,col));
                
                if (matrix_elemet > max_val)
                {   // saves value and position if it is larger than the
                    // previous max value
                    max_val = matrix_elemet;
                    idx_row = row;
                    idx_col = col;
                }
            }
        }
    }
    
    return max_val;
}


void transform(int n, arma::mat& A, arma::mat& R, int idx_row, int idx_col)
{
    /*
    This function rotatetes the (symmetric) matrix A by an angle theta in the 
    plane spanned by the unit vectors e_l and e_k, where 
    e_i = [0, 0, ..., 0, 1, 0, ..., 0]. The angle theta is defined s.t. the 
    elements a_lk = a_kl = 0 after the transformation.

    Parameters
    ----------
    A : arma::mat&
        Reference to matrix A.

    R : arma::mat&
        Reference to eigenvector matrix R.

    idx_row : int
        The index k of one of the unit vectors that are used in the rotation.

    idx_col : int
        The index l of the other unit vector that are used in the rotation.
    */



    double a_kk = A(idx_row, idx_row);
    double a_ll = A(idx_col, idx_col);
    double a_lk = A(idx_col, idx_row);
    double r_ik;
    double r_il;

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
    
    double c = 1/std::sqrt(1 + t*t);
    double s = t*c;
    double a_ki;

    for (int i = 0; i < n; i++)
    {   // elements kk, ll, lk, kl should not be changed, but it is fixed after
        // the loop instead of using an if-test in the loop (vectorization)
        
        a_ki = A(idx_row, i);    // since this matrix element will be overwritten
        
        A(i, idx_row) = A(idx_row, i) = a_ki*c - A(idx_col, i)*s;
        A(i, idx_col) = A(idx_col, i) = A(idx_col, i)*c + a_ki*s;

        r_ik = R(i, idx_row);
        r_il = R(i, idx_col);

        R(i, idx_row) = c*r_ik - s*r_il;
        R(i, idx_row) = c*r_il + s*r_ik;
    }

    A(idx_col, idx_row) = A(idx_row, idx_col) = 0; // by virtue of the algorithm (this is how we found theta!).
    A(idx_row, idx_row) = a_kk*c*c - 2*a_lk*c*s + a_ll*s*s;
    A(idx_col, idx_col) = a_ll*c*c + 2*a_lk*c*s + a_kk*s*s;

}


void find_eig(int n, arma::mat& A, double tol_off_diag)
{
    /*
    Finds the eigenvalues of matrix A. Uses find_max to locate the largest value
    in the array. Uses transform to eliminate off-diagonal elements. Repeats the
    process untill all diagonal elements are smaller than tol_off_diag.

    Creates a matrix R, for storing eigenvectors.

    Parameters
    ----------
    n : int
        Dimension of array.

    A : arma::mat&
        Reference to matrix.

    tol_off_diag : double
        Tolerance for the largest allowed value of the off-diagonal elements.
    */

    arma::mat R(n, n);      // matrix for storing the eigenvectors of A
    R.zeros();
    R.diag(0) += 1;
    
    int idx_col;
    int idx_row;
    double max_val = find_max(n, A, idx_col, idx_row);
    
    while (max_val > tol_off_diag)
    {
        transform(n, A, R, idx_col, idx_row);
        max_val = find_max(n, A, idx_col, idx_row);
    }

    R.print();
}


void test_find_max()
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


    // std::cout << max_val << std::endl;  // debug, remove later
}


void test_inner_product_conserved()
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
    A.print();

    for (int i = 0; i < n; i++)
    {   // loops over all column vectors
        
        for (int j = 0; j < n; j++)
        {   // loops over all column vectors
            std::cout << arma::dot(A.col(i), A.col(j)) << std::endl;
        }
    }
}


void test_find_eig()
{
    /*
    Checks that the find_eig function finds eigenvalues up to a fixed tolerance.
    */
    
    double pi = 3.14159265358979323846;
    int n = 5;                       // dimension of matrix
    double step     = 1.0/n;
    double diag     = 2/(step*step);    // diagonal elements
    double off_diag = -1/(step*step);   // off-diagonal elements
    double tol_off_diag = std::pow(10, -15);    // tolerance of off-diagonal elements
    double tol_eig  = 1;                // tolerance of eigenvalues

    // std::cout << diag << std::endl;
    
    arma::mat A = construct_diag_matrix(n, off_diag, diag);
    // arma::mat A = construct_diag_matrix(n);
    arma::vec eigenvalues(n);

    
    for (int i = 1; i <= n; i++)
    {   // generating analytical eigenvalues
        eigenvalues(i - 1) = diag + 2*off_diag*std::cos(i*pi/(n + 1));
    }
    
    find_eig(n, A, tol_off_diag); // numerical eigenvalues

    arma::vec sorted_diag = arma::sort(A.diag(0));
    
    for (int i = 0; i < n; i++)
    {   // checking that analytical and numerical results match to a given tolerance
        
        if (fabs(eigenvalues(i) - sorted_diag(i)) > tol_eig)
        {   
            std::cout << "ERROR: Exact eigenvalue = " 
                << eigenvalues(i) << ", computed eigenvalue = "
                << sorted_diag(i) << "." << std::endl;
        }

    }
    // A.print();
}


int main(int argc, char* argv[]) {
    // int n = 5;
    // double stepsize = 1/(n+1);

    // arma::mat A = construct_diag_matrix(n);
    // A.print();
    test_find_max();
    // test_inner_product_conserved();
    test_find_eig();


    return 1;
}