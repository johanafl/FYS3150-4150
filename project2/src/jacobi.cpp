#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <chrono>


arma::mat construct_diag_matrix(int n) {
    /*
    Function for creating a diagonal armadillo matrix.
    */
    arma::mat A(n, n);
    A.zeros();

    A.diag(0) += 2;
    A.diag(-1) += -1.0;
    A.diag(1) += -1.0;

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

double find_max(int n, arma::mat A, int* idx_row, int* idx_column) {
    /*
    Finds the max value of the off-diagonal elements of the matrix. Saves
    the indexes of the max value.

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
            {
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

void transform(int n, arma::mat* A, int idx_col, int idx_row) {
    /*
    This function rotatetes the (symmetric) matrix A by an angle theta in the 
    plane spanned by the unit vectors e_l and e_k, where 
    e_i = [0,0,...,0,1,0,...,0]. The angle theta is defined s.t. that the 
    elements a_lk = a_kl = 0 after the transformation.

    Parameters
    ----------
    A : arma::mat*
        Pointer to the symmetric matrix A
    idx_row : int
        The index k of one of the unit vectors that are used in the rotation.
    idx_col : int
        The index l of the other unit vector that are used in the rotation.
    */
    arma::mat B = *A;
    double a_ll = B(idx_row, idx_row);
    double a_kk = B(idx_col, idx_col);
    double a_kl = B(idx_col, idx_row);

    double tau = (a_ll - a_kk)/(2*a_kl); // We have defined tau = cot(2*theta), where theta is unknown. We are choosing thata s.t. the element b_kl = b_lk = 0, where B is the new matrix after the transformation with elements b_ij.
    double t = -tau + std::sqrt(1 + tau*tau); // We have defined t = tan(theta). We could have choosen minus instead of pluss in front of the square root, and gotten the same result (by virtue of the algorithm).
    double c = 1/std::sqrt(1 + t*t);
    double s = t*c;

    for (int i=0; i<n; i++) // NB! In project description: idx_row = k and idx_col = l.
    {
        B(i,idx_row) = B(idx_row,i) = B(idx_row,i)*c - B(idx_col,i)*s;
        B(i,idx_col) = B(idx_col,i) = B(idx_col,i)*c + B(idx_row,i)*s;
    }    
    // B(idx_col,idx_row) = B(idx_row,idx_col) = 0; // By virtue of the algorithm (this is how we found theta!).
    B(idx_row,idx_row) = B(idx_row,idx_row)*c*c - 2*B(idx_col,idx_row)*c*s + B(idx_col,idx_col)*s*s;
    B(idx_col,idx_col) = B(idx_col,idx_col)*c*c - 2*B(idx_col,idx_row)*c*s + B(idx_row,idx_row)*s*s;
}

arma::mat find_eig(int n, arma::mat A, double eps) {
    /*
    Not implemented yet.
    */
    arma::mat B = A;
    int idx_col;
    int idx_row;
    double max_val = find_max(n, B, & idx_col, & idx_row);
    while (max_val > eps)
    {
        transform(n, & B, idx_col, idx_row);
        max_val = find_max(n, B, & idx_col, & idx_row);
    }
    return B;
}

void test_find_max() {
    /*
    Test function that tests the most important functionality of the Jacobi class.
    */
    int n = 5;
    double maximum = 5555;
    arma::mat A(n,n);
    A.zeros();
    for (int i=0; i<5; i++)
    {
        for (int j=0; j<5; j++)
        {
            A(i,j) = i*j*1.5;
        }
    }
    A(3,2) = maximum;
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
}


void test_inner_product_conserved() {
    /*
    Checks that the inner product is preserved after one rotation with the
    transformation matrix.
    */
    int n = 3;
    arma::mat A(n,n);
    A.zeros();
    A(0, 2) = 1;
    A(1, 0) = 1;
    A(2, 1) = 1;
    
    double tol = 1;
    int idx_row = 0;
    int idx_column = 2;

    transform(n, & A, idx_row, idx_column); // This must be implemented.
    double inner_prod = arma::dot(A.col(0),A.col(1));
    
    if (fabs(inner_prod) > tol)
    {
        std::cout << "Something went wrong!" << std::endl;
    }
}

void test_find_eig() {
    /*
    Checks that the find_eig function finds eigenvalues up to a set error eps.
    */
    double eps_eig = 0.01;
    double eps_test = 1;
    int n = 5;
    arma::mat A = construct_diag_matrix(n);
    arma::vec lam(n);
    double d = 2;
    double a = -1;
    double pi = 3.14159265358979323846;
    for (int i=0; i<n; i++)
    {
        lam(i) = d + 2*a*std::cos(i*pi/(n+1));
    }
    
    arma::mat B = find_eig(n, A, eps_eig);
    for (int i=0; i<n; i++)
    {
        if (fabs(lam(i) - B(i,i)) > eps_test)
        {
            std::cout << "Something went wrong!\n" << "Exact eigenvalue = " 
                      << lam(i) << ", computed eigenvalue = " << B(i,i)
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
    // test_inner_product_conserved();
    test_find_eig();

    return 1;
}