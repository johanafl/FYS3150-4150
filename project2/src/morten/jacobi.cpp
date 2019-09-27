#include "jacobi.h"

void jac::hei() {
    std::cout << "hei" << std::endl;
}

arma::mat jac::construct_diag_matrix(int n)
{
    /*
    Function for creating a tridiagonal matrix. Second derivative discretized.

    CURRENTLY NOT USED

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


arma::mat jac::construct_diag_matrix(int n, double off_diag, double diag)
{
    /*
    Function for creating a tridiagonal matrix. With inputs for
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

arma::mat jac::construct_diag_matrix(int n, double off_diag, arma::vec diag)
{
    /*
    Function for creating a tridiagonal matrix. With inputs for
    populating the upper and lower diagonal where diagonal input is a vector.

    Parameters
    ----------
    n : int
        Dimension of matrix.

    off_diag : double
        Value for the upper and lower diagonal elements.

    diag : arma::vec
        Values for the diagonal elements.
    */

    arma::mat A(n, n);
    
    A.zeros();
    A.diag(0)  = diag;
    A.diag(-1) += off_diag;
    A.diag(1)  += off_diag;

    return A;
}


double jac::find_max(int n, arma::mat& A, int& idx_row, int& idx_col)
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


void jac::transform(int n, arma::mat& A, arma::mat& R, int idx_row, int idx_col)
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


void jac::find_eig(int n, arma::mat& A, double tol_off_diag)
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
}