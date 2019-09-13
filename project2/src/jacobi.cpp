#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <chrono>


class Jacobi {
/*
Class for finding eigenvalues of matrix using Jacobi's method.
*/
private:
    // arma::mat orto_transform_mat;
    arma::mat mat;
    int n;
    double eps;

public:
    Jacobi(int size, arma::mat original_mat, double tolorance) 
    {
    /*
    Constructer of the Jacobi class. 
    */
        n   = size;
        eps = tolorance;
        mat = original_mat;

        arma::mat orto_transform_mat(n,n);
        orto_transform_mat.diag(0) += 1.0;
    }

    void create_transformation_matrix(double theta, int l)
    {
    /*
    Constructer of the Jacobi class. 
    */
        double c = std::cos(theta);
        double s = std::sin(theta);
        orto_tranformation_mat(l,l) = c;
        orto_tranformation_mat(l,n) = s;
        orto_tranformation_mat(n,l) = -s;
        orto_tranformation_mat(n,n) = c;
    }

    void step_forward()
    {
    /*
    Perhaps call the create_transformation_matrix
    */
        // Compute next step
    }

    void find_max()
    {
        double max_val = 0;
        int idx_l;
        int idx_k;

        for (int i=0; i<n; i++)
        {
            for (int j=0; j<n; j++)
            {
                if (mat(i,j) > max_val)
                {
                    max_val = 
                }
            }
        }
    }

    void find_eigenvalues ()
    {
    /*
    Check if norm is less than epsilon. If not, continue to find eigenvalues
    */
        // implement a while loop. norm(mat,"fro") should take the frobenius norm (http://arma.sourceforge.net/docs.html)
        while (norm(mat,"fro") > eps)
        {
            step_forward()
        }
    }

    ~Jacobi() 
    {
    /*
    Deallocating all allocated objects.
    */
    // delete[] orto_transform_mat;
    }
};


void test() {
  /*
  Test function that tests the most important functionality of the Jacobi class.
  */
}


int main(int argc, char* argv[])
{
    return 1;
}