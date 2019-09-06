#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <chrono>

void write_to_file(std::string filename, int n, double*v);
double relative_error(double v, double u);
void LU_arma(int n);

double rhs_func(double x)
/*
Function for calculating the right-hand-side of the equation, i.e.,
f(x) = 100e^(-10x).
*/
{
    return 100*exp(-10*x);
}

double exact_solution(double x)
/*
Function for calculating the exact solution to our numerical problem. This will
be used for comparing the numerical precision for our computed solution.
*/
{
    return 1 - (1 - exp(-10))*x - exp(-10*x);
}


void thomas_algorithm(int n) // include diaginal variables
/*
Function for doing Gaussian elimination on our trigonal matrix, i.e.,
Thomas algorithm, and saving the the values to a .txt-file.

Since we have a trigonal matrix, this function will iterate over three vectors
instead of the whole nxn matrix; one vector with n elements which will
represent the diagonal elements, and two vectors with n-1 elements which will
represent the upper and lower diagonal elements.
This function will print out the time used by the gaussian elimination
algorithm.
*/
{
    double stepsize = 1.0/(n+1);

    double* lower_diag = new double[n-1]; 
    double* diag       = new double[n];   
    double* upper_diag = new double[n-1]; 
    double* rhs_val    = new double[n];   
    double* computed   = new double[n];   

    for (int i=0; i<n; i++)
    {
        double x = stepsize*(i+1);
        rhs_val[i] = rhs_func(x)*(stepsize*stepsize);
    }

    for (int i=0; i<n-1; i++)
    {
        lower_diag[i] = -1.0;
        diag[i]       = 2.0;
        upper_diag[i] = -1.0;
    }
    diag[n-1] = 2.0;

    ////////////////////////////////
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now(); // Taken from this site: http://www.cplusplus.com/reference/chrono/steady_clock/
    ////////////////////////////////

    for (int i=0; i<n-1; i++)
    {
        diag[i+1] = diag[i+1] - lower_diag[i]/diag[i]*upper_diag[i];
        rhs_val[i+1] = rhs_val[i+1] - lower_diag[i]/diag[i]*rhs_val[i];
    }

    computed[n-1] = rhs_val[n-1]/diag[n-1];
    for (int i=n-1; i>=1; i--)
    {
        computed[i-1] = (rhs_val[i-1] - upper_diag[i-1]*computed[i])/diag[i-1];
    }

    ////////////////////////////////
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> tot_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    ////////////////////////////////
    
    std::cout << std::setprecision(32) << tot_time.count() << std::endl; // NB! Gives time in seconds.

    std::string filename = "Poisson_values_n_" + std::to_string(n) + ".txt";
    write_to_file(filename, n, computed);

    delete[] lower_diag;
    delete[] diag;
    delete[] upper_diag;
    delete[] rhs_val;
    delete[] computed;
}

void thomas_algorithm_special(int n)
/*
Function for doing Gaussian elimination for the special trigonal matrix, i.e.,
Thomas algorithm on a special trigonal matrix, and saving the the values to a
.txt-file.

This function works just like the general Thomas algorithm function, only here
we don't need to include the upper and lower diagonal elements of the matrix
since they er both just one, saving us some FLOPS.
*/
{
    double stepsize = 1.0/(n+1);

    double* diag     = new double[n];   
    double* rhs_val  = new double[n];   
    double* computed = new double[n];   

   for (int i=0; i<n; i++)
    {
        double x = stepsize*(i+1);
        rhs_val[i] = rhs_func(x)*(stepsize*stepsize);
    }

    for (int i=0; i<n; i++)
    {
        diag[i] = 2.0;
    }

    ////////////////////////////////
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now(); // Taken from this site: http://www.cplusplus.com/reference/chrono/steady_clock/
    ////////////////////////////////

    for (int i=0; i<n-1; i++)
    {
        diag[i+1] = 2.0 - 1.0/diag[i];
        rhs_val[i+1] = rhs_val[i+1] + 1.0/diag[i]*rhs_val[i];
    }
    // Counting 2 + 3 = 5 FLOPS. Is this correct? (Why 4?)

    computed[n-1] = rhs_val[n-1]/diag[n-1];
    for (int i=n-1; i>=1; i--)
    {
        computed[i-1] = (rhs_val[i-1] + computed[i])/diag[i-1];
    }
    // Counting 2 + 3 = 5 FLOPS. Is this correct?

    ////////////////////////////////
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> tot_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    ////////////////////////////////
    
    std::cout << std::setprecision(32) << tot_time.count() << std::endl; // NB! Gives time in seconds.

    std::string filename = "Poisson_values_spes_n_" + std::to_string(n) + ".txt";
    write_to_file(filename, n, computed);

    delete[] diag;
    delete[] rhs_val;
    delete[] computed;
}

void write_to_file(std::string filename, int n, double*computed_val)
/*
Function for writing the values to a .txt-file.
Writing the exact and computed solution as well as the error for the computed
relative to the exact solution.
*/
{
    double* exact_val = new double[n];
    double* eps       = new double[n];

    double stepsize    = 1.0/(n+1);

    for (int i=0; i<n; i++)
    {
        double x = stepsize*(i+1);
        exact_val[i]     = exact_solution(x);
        if(exact_val[i] == 0)
        {
            eps[i] = -1;
        }
        else
        {
            eps[i] = relative_error(computed_val[i], exact_val[i]);
        }    
    }

    std::ofstream results_file;
    results_file.open(filename);
    results_file << std::setw(25) << "U(x) (exact)"
                 << std::setw(25) << "V(x) (discretized)"
                 << std::setw(25) << "Relative error \n"
                 << std::setw(25) << 0
                 << std::setw(25) << 0
                 << std::setw(25) << 0 << std::endl;
    for (int i=0; i<n; i++)
    {
        results_file << std::setw(25) << std::setprecision(16) << exact_val[i]
                     << std::setw(25) << std::setprecision(16) << computed_val[i]
                     << std::setw(25) << std::setprecision(16) << eps[i]
                     << std::endl;
    }
    results_file << std::setw(25) << 0
                 << std::setw(25) << 0
                 << std::setw(25) << 0 << std::endl;
    results_file.close();
}

void write_to_file(std::string filename, int n, arma::vec computed_val)
/*
Function for writing the values to a .txt-file.
Writing the exact and computed solution as well as the error for the computed
relative to the exact solution.
This function is used when we use armadillo to do the Thomas algorithm.
*/
{
    double* computed = new double[n];
    for(int i=0; i<n; i++)
    {
        computed[i] = computed_val[i];
    }
    write_to_file(filename, n, computed);
    delete[] computed;
}

double relative_error(double computed_val, double exact_val)
/*
Function for calculating the relative error for the computed solution relative
to the exact solution.
*/
{
    return fabs((computed_val - exact_val)/exact_val);
}

void LU_arma(int n)
/*
Function for calculating the Thomas algorithm using the armadillo library.

This function also prints the time it takes to run the algorithm.
*/
{
    double stepsize = 1.0/(n+1);

    arma::mat A(n, n);

    A.diag(0)  += 2.0;
    A.diag(1)  += -1.0;
    A.diag(-1) += -1.0;

    arma::vec rhs_val(n);
    arma::vec tmp(n);
    arma::vec computed;

    arma::mat L;
    arma::mat U;

   for (int i=0; i<n; i++)
    {
        double x = stepsize*(i+1);
        rhs_val(i) = rhs_func(x)*(stepsize*stepsize);
    }

    ////////////////////////////////
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now(); // Taken from this site: http://www.cplusplus.com/reference/chrono/steady_clock/
    ////////////////////////////////

    arma::lu(L, U, A);

    arma::solve(tmp, L, rhs_val);
    arma::solve(computed, U, tmp);

    ////////////////////////////////
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> tot_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    ////////////////////////////////
    std::cout << std::setprecision(32) << tot_time.count() << std::endl; // NB! Gives time in seconds.

    std::string filename = "Poisson_values_LU_n_" + std::to_string(n) + ".txt";
    write_to_file(filename, n, computed);
}

int main(int argc, char *argv[]) 
{
    int n = atoi(argv[1]);
    thomas_algorithm(n);
    thomas_algorithm_special(n);
    LU_arma(n);
 
    return 1;
}