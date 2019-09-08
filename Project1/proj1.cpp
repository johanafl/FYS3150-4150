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

Parameters
----------
x : double
    Function variable.
*/
{
    return 1 - (1 - exp(-10))*x - exp(-10*x);
}


double thomas_algorithm(int n, bool write)
/*
Function for doing Gaussian elimination on our trigonal matrix, i.e.,
Thomas algorithm, and saving the the values to a .txt-file.

Since we have a trigonal matrix, this function will iterate over three vectors
instead of the whole nxn matrix; one vector with n elements which will
represent the diagonal elements, and two vectors with n-1 elements which will
represent the upper and lower diagonal elements.

Parameters
----------
n : int
    Dimension of matrix. Number of discrete points.

write : bool
    Boolean for toggling write to file on/off.
*/
{
    double stepsize = 1.0/(n+1);

    double* lower_diag = new double[n-1]; 
    double* diag       = new double[n];   
    double* upper_diag = new double[n-1]; 
    double* rhs_val    = new double[n];   
    double* computed   = new double[n];   

    for (int i=0; i<n; i++)
    {   // calculating the r.h.s. values
        double x = stepsize*(i+1);
        rhs_val[i] = rhs_func(x)*(stepsize*stepsize);
    }

    for (int i=0; i<n-1; i++)
    {   // inserting values in the diagonals
        lower_diag[i] = -1.0;
        diag[i]       = 2.0;
        upper_diag[i] = -1.0;
    }
    diag[n-1] = 2.0;

    // Copied from: http://www.cplusplus.com/reference/chrono/steady_clock/
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

    for (int i=0; i<n-1; i++)
    {   // elimination of bottom diagonal by forward substitution
        diag[i+1] = diag[i+1] - lower_diag[i]/diag[i]*upper_diag[i];
        rhs_val[i+1] = rhs_val[i+1] - lower_diag[i]/diag[i]*rhs_val[i];
    }

    computed[n-1] = rhs_val[n-1]/diag[n-1];
    for (int i=n-1; i>=1; i--)
    {   // elimination of top diagonal by backwards substitution
        computed[i-1] = (rhs_val[i-1] - upper_diag[i-1]*computed[i])/diag[i-1];
    }

    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> total_time = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
    

    if (write) 
    {
        std::cout << std::setprecision(32) << total_time.count() << std::endl; // NB! Gives time in seconds.
        std::string filename = "Poisson_values_n_" + std::to_string(n) + ".txt";
        write_to_file(filename, n, computed);
    }

    return computed;

    delete[] lower_diag;
    delete[] diag;
    delete[] upper_diag;
    delete[] rhs_val;
    delete[] computed;

    return total_time.count();
}

double thomas_algorithm_special(int n, bool write)
/*
Function for doing Gaussian elimination for the special trigonal matrix, i.e.,
Thomas algorithm on a special trigonal matrix, and saving the the values to a
.txt-file.

This function works just like the general Thomas algorithm function, only here
we don't need to include the upper and lower diagonal elements of the matrix
since they er both just one, saving us some FLOPS.

Parameters
----------
n : int
    Dimension of matrix. Number of discrete points.
*/
{
    double stepsize = 1.0/(n+1);

    double* diag     = new double[n];   
    double* rhs_val  = new double[n];   
    double* computed = new double[n];   

   for (int i=0; i<n; i++)
    {   // calculating the r.h.s.
        double x = stepsize*(i+1);
        rhs_val[i] = rhs_func(x)*(stepsize*stepsize);
    }

    for (int i=0; i<n; i++)
    {   // inserting values in the diagonal
        diag[i] = 2.0;
    }

    for (int i=0; i<n-1; i++)
    {   // forward substituting, separate loop to exclude from timing
        diag[i+1] = 2.0 - 1.0/diag[i];
    }

    // starting timer
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

    for (int i=0; i<n-1; i++)
    {   // forward substituting
        rhs_val[i+1] = rhs_val[i+1] + rhs_val[i]/diag[i];
    }

    computed[n-1] = rhs_val[n-1]/diag[n-1];
    for (int i=n-1; i>=1; i--)
    {   // backward substituting
        computed[i-1] = (rhs_val[i-1] + computed[i])/diag[i-1];
    }

    // ending timer
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> total_time = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);



    if (write) 
    {   
        std::cout << std::setprecision(32) << total_time.count() << std::endl; // NB! Gives time in seconds.
        std::string filename = "Poisson_values_spes_n_" + std::to_string(n) + ".txt";
        write_to_file(filename, n, computed);
    }

    return computed;

    delete[] diag;
    delete[] rhs_val;
    delete[] computed;

    return total_time.count();
}

void write_to_file(std::string filename, int n, double* computed_val)
/*
Function for writing the values to a .txt-file.
Writing the exact and computed solution as well as the error for the computed
relative to the exact solution.

Parameters
----------
filename : std::string
    A string input with the name of the file which data are written to.

n : int
    Dimension of matrix. Number of discrete points.
*/
{
    
    double* exact_val = new double[n];
    double* eps       = new double[n];
    double stepsize   = 1.0/(n+1);
    double x;

    for (int i=0; i<n; i++)
    {   // calculating exact solution and relative error
        x = stepsize*(i+1);
        exact_val[i] = exact_solution(x);
        eps[i] = relative_error(computed_val[i], exact_val[i]);
       
    }

    // writing titles and first known values to file
    std::ofstream results_file;
    results_file.open(filename);
    results_file << std::setw(25) << "U(x) (exact)"
                 << std::setw(25) << "V(x) (discretized)"
                 << std::setw(25) << "Relative error \n"
                 << std::setw(25) << 0
                 << std::setw(25) << 0
                 << std::setw(25) << 0 << std::endl;
    for (int i=0; i<n; i++)
    {   // writing calculated values to file
        results_file << std::setw(25) << std::setprecision(16) << exact_val[i]
                     << std::setw(25) << std::setprecision(16) << computed_val[i]
                     << std::setw(25) << std::setprecision(16) << eps[i]
                     << std::endl;
    }

    // writing final known values to file
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

void write_to_file_error(std::string filename, int n, double*computed_val)
/*
Function for writing the relative error to a .txt-file.
*/
{
    double exact_val;
    double eps = 0;

    double stepsize    = 1.0/(n+1);

    for (int i=0; i<n; i++)
    {
        double x = stepsize*(i+1);
        exact_val     = exact_solution(x);
        if(eps < relative_error(computed_val[i], exact_val))
        {
            eps = relative_error(computed_val[i], exact_val);
        }
    }

    std::ofstream error_file;
    error_file.open(filename, std::ios_base::app);
    error_file << std::setw(25) << std::setprecision(16) << eps << std::endl;
    error_file.close();
}

double relative_error(double computed_val, double exact_val)
/*
Function for calculating the relative error for the computed solution relative
to the exact solution.
*/
{
    return fabs((computed_val - exact_val)/exact_val);
}

double LU_arma(int n, bool write)
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

    L.zeros();
    U.zeros();

   for (int i=0; i<n; i++)
    {   // calculating r.h.s. values
        double x = stepsize*(i+1);
        rhs_val(i) = rhs_func(x)*(stepsize*stepsize);
    }

    // starting timer
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

    arma::lu(L, U, A);
    arma::solve(tmp, L, rhs_val);
    arma::solve(computed, U, tmp);

    // ending timer
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> total_time = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);


    if (write)
    {
        std::cout << std::setprecision(32) << total_time.count() << std::endl; // NB! Gives time in seconds.
        std::string filename = "Poisson_values_LU_n_" + std::to_string(n) + ".txt";
        write_to_file(filename, n, computed);
    }

    return total_time.count();


}

void compare_times()
/*
Function for comparing computation times of thomas, thomas special, and LU.
Each grid size is calculated 10 times, and all timings are written to a text
file compare_times.txt.
*/
{
    int grid_values = 16;    // number of different grid values
    int N[grid_values];
    int runs = 10;          // number of runs for each grid size

    N[0]  = 10;
    N[1]  = 100;
    N[2]  = 500;
    N[3]  = 1000;
    N[4]  = 5000;
    N[5]  = 10000;
    N[6]  = 50000;
    N[7]  = 100000;
    N[8]  = 500000;
    N[9]  = 1000000;
    N[10] = 5000000;
    N[11] = 10000000;
    N[12] = 50000000;
    N[13] = 100000000;
    N[14] = 500000000;
    N[15] = 1000000000;


    // creating file and writing headers
    std::string filename = "compare_times.txt";
    std::ofstream compare_times_file;
    compare_times_file.open(filename);
    compare_times_file  << std::setw(40) << "thomas algorithm"
                        << std::setw(40) << "thomas algorithm special"
                        << std::setw(40) << "LU\n"
                        << std::setw(40) << runs
                        << std::setw(40) << runs
                        << std::setw(40) << runs << "\n"
                        << std::setw(40) << grid_values
                        << std::setw(40) << grid_values
                        << std::setw(40) << grid_values << "\n";

    for (int i=0; i<grid_values; i++)
    {   // looping over each grid size
        // writing grid size info to file
        std::cout << "calculating grid size " + std::to_string(N[i])
            + " of " + std::to_string(N[grid_values-1]) << std::endl;
        compare_times_file  << std::setw(40) << std::to_string(N[i])
                            << std::setw(40) << std::to_string(N[i])
                            << std::setw(40) << std::to_string(N[i])
                            << "\n";
        
        for (int _=0; _<runs; _++)
        {   // looping over each grid size 'runs' amount of times
            // writing timing data to file
            compare_times_file  << std::setw(40) << std::setprecision(32) 
                                << thomas_algorithm(N[i], false)
                                << std::setw(40) << std::setprecision(32)
                                << thomas_algorithm_special(N[i], false);
            
            if (N[i] <= 10000)
            {   // limits the gid size for LU since the calculations are
                // impossible with normal hardware at values approaching 100000
            compare_times_file  << std::setw(40) << std::setprecision(32)
                                << LU_arma(N[i], false) << std::endl;
            }
            
            else
            {   // writes -1 if LU is incapable of providing results
            compare_times_file  << std::setw(40) << std::setprecision(32)
                                << -1 << std::endl;  
            }
                                
        }
    }

    compare_times_file.close();
}


int main(int argc, char *argv[])
{
    int n = atoi(argv[1]);
    // thomas_algorithm(n);
    // thomas_algorithm_special(n);
    // LU_arma(n);
    // compare_times();
 
    return 1;
}