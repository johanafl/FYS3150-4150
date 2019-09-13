#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <chrono>



void create_transformation_matrix(double theta, int l) {
/*
Not implemented yet. Supposed to create an orthogonal transformation matrix.
Think this function is unnessesary.
*/
}

void step_forward() {
/*
Not implemented yet. Perhaps call the create_transformation_matrix. Supposed to 
multiply the transformation matrix with the matrix we have. 
*/
    // Compute next step
}


void find_eigenvalues() {
/*
Not implemented yet. Check if norm is less than epsilon. If not, continue to 
find eigenvalues
*/
    // (old)    implement a while loop. norm(mat,"fro") should take the frobenius norm (http://arma.sourceforge.net/docs.html)
    // (update) We are supposed to use the maximum value. use the find_max function.

    // while ( > tolerance)
    // {
    //     step_forward();
    // }
}



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

    return A;
}

double find_max(int n, arma::mat A, int idx_row, int idx_column) {
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
                    idx_row = row;
                    idx_column = column;
                }
            }
        }
    }
    return max_val;
}

void test() {
    /*
    Test function that tests the most important functionality of the Jacobi class.
    */

}



int main(int argc, char* argv[]) {
    int n = 5;
    double stepsize = 1/(n+1);

    arma::mat A = construct_diag_matrix(n);
    A.print();
    return 1;
}