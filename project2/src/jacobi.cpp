#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <chrono>


// class Jacobi {
// /*
// Class for finding eigenvalues of matrix using Jacobi's method.
// */
// private:
//     // arma::mat orto_transform_mat;
//     arma::mat mat;
//     int n;
//     double eps;

// public:
//     Jacobi(int size, arma::mat original_matrix, double tolerance) {
//     /*
//     Constructer of the Jacobi class. 
//     */

//         arma::mat orto_transform_mat(size, size);
//         orto_transform_mat.zeros();
//         orto_transform_mat.diag(0) += 1.0;
//     }

//     void create_transformation_matrix(double theta, int l) {
//     /*

//     */
    
//     }

//     void step_forward() {
//     /*
//     Perhaps call the create_transformation_matrix
//     */
//         // Compute next step
//     }

    // void find_max() {
    //     /*
    //     Finds the max value of the off-diagonal elements of the matrix. Saves
    //     the indexes of the max value.
    //     */

    //     double max_val = 0;     // final max value
    //     double check_max;       // placeholder for new max value to check
    //     int max_row;            // row position of max value
    //     int max_column;         // column position of max value

    //     for (int row=0; i<n; i++)
    //     {   // loops over rows
            
    //         for (int column=0; j<n; j++)
    //         {   // loops over columns
    //             check_max = original_matrix(i,j);
                
    //             if (check_max > max_val)
    //             {
    //                 max_val = check_max;
    //                 max_row = row;
    //                 max_column = column;
    //             }
    //         }
    //     }
    // }

//     void find_eigenvalues() {
//     /*
//     Check if norm is less than epsilon. If not, continue to find eigenvalues
//     */
//         // implement a while loop. norm(mat,"fro") should take the frobenius norm (http://arma.sourceforge.net/docs.html)
//         while (norm(mat,"fro") > tolerance)
//         {
//             step_forward();
//         }
//     }

//     ~Jacobi() {
//     /*
//     Deallocating all allocated objects.
//     */
//     // delete[] orto_transform_mat;
//     }
// };


void test() {
    /*
    Test function that tests the most important functionality of the Jacobi class.
    */

}

void construct_matrix() {

    int step;
    int n = 10;
    arma::mat A(n, n);
    A.zeros();

    A.diag(0) += 1;
    A.diag(-1) += -1.0/(step*step);

}


int main(int argc, char* argv[]) {

    construct_matrix();
    return 1;
}