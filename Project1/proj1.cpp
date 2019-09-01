#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include "time.h"   // Things from this does not work. Ask next time!
#include <chrono>

void write_to_file(std::string filename, int n, double*v);
double relative_error(double v, double u);

double rhs_func(double x)
{
    return 100*exp(-10*x);
}

double exact_solution(double x)
{
    return 1 - (1 - exp(-10))*x - exp(-10*x);
}


void gauss_elim(int n) // include diaginal variables
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
        // std::cout << x << std::endl;
        rhs_val[i] = rhs_func(x)*(stepsize*stepsize);
    }
    // std::cout << std::endl;
    // for (int i=0; i<n; i++)
    // {
    //     std::cout << rhs_val[i] << std::endl;
    // }

    for (int i=0; i<n-1; i++)
    {
        lower_diag[i] = -1.0;
        diag[i]       = 2.0;
        upper_diag[i] = -1.0;
    }
    diag[n-1] = 2.0;

    ////////////////////////////////
    // clock_t start, finish; // Gives only integer time. Not double!
    // start = clock();

    // time_t time0; // create timers.   // This also give wrong.
    // time_t time1;
    // time(&time0); // get current time.

    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now(); // Taken from this site: http://www.cplusplus.com/reference/chrono/steady_clock/
    ////////////////////////////////

    for (int i=0; i<n-1; i++)
    {
        // std::cout << lower_diag[i] << std::endl;
        // std::cout << upper_diag[i] << std::endl;
        diag[i+1] = diag[i+1] - lower_diag[i]/diag[i]*upper_diag[i];
        rhs_val[i+1] = rhs_val[i+1] - lower_diag[i]/diag[i]*rhs_val[i];
        // std::cout << diag[i+1] << std::endl;
    }

    computed[n-1] = rhs_val[n-1]/diag[n-1];
    for (int i=n-1; i>=1; i--)
    {
        // std::cout << lower_diag[i-1] << std::endl;
        // std::cout << upper_diag[i-1] << std::endl;
        // std::cout << diag[i-1] << std::endl;
        // computed[i-1] = (rhs_val[i-1] - upper_diag[i-1]*computed[i])/diag[i-1];
        computed[i-1] = (rhs_val[i-1] - upper_diag[i-1]*computed[i])/diag[i-1];
    }
    ////////////////////////////////
    // finish = clock();
    // double tot_time = (finish - start)/CLOCKS_PER_SEC;

    // time(&time1);   // get current time after time pass.
    // double tot_time = time1 - time0;
 
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

/*
// Old solution:
void gauss_elim(int n) 
{
    double h = 1.0/(n+1);

    double* a      = new double[n-1]; 
    double* b      = new double[n];   
    double* c      = new double[n-1]; 
    double* f      = new double[n];   
    double* v      = new double[n];   
    double* f_prim = new double[n];   
    double* b_prim = new double[n];   

   for (int i=0; i<n; i++)
    {
        double x = h*(i+1);
        f[i] = func_f(x);
    }


    for (int i=0; i<n-1; i++)
    {
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;
    }
    b[n-1] = 2;

    b_prim[0] = b[0];
    for (int i=0; i<n-1; i++)
    {
        b_prim[i+1] = b[i+1] - a[i]/b[i]*c[i];
        f_prim[i+1] = f[i+1] - a[i]/b[i]*f[i];
    }
    // v[n-1] = f_prim[n-1];
    // for (int i=n-1; i>=1; i--)
    // {
    //     v[i-1] = f_prim[i-1] - c[i-1]*v[i];
    // }
    for (int i=n-1; i>=1; i--)
    {
        f_prim[i-1] = f_prim[i-1] - c[i-1]/b_prim[i]*f_prim[i];
    }
    for (int i=0; i<n; i++)
    {
        v[i] = f_prim[i]/b_prim[i];
    }
    
    std::string filename = "Poisson_values_n_" + std::to_string(n)
                             + ".txt";
    write_to_file(filename, n, v);

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] f;
    delete[] v;
    delete[] f_prim;
    delete[] b_prim;
}
*/

void gauss_elim_spes(int n) 
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

    // Maybe change the name of the file to indecate that the specialized function was used?
    std::string filename = "Poisson_values_spes_n_" + std::to_string(n) + ".txt";
    write_to_file(filename, n, computed);

    delete[] diag;
    delete[] rhs_val;
    delete[] computed;
}

/*
void gauss_elim_spes(int n) 
{
    double h = 1.0/(n+1);

    double* a      = new double[n-1]; 
    double* b      = new double[n];   
    double* c      = new double[n-1]; 
    double* f      = new double[n];   
    double* v      = new double[n];   
    double* f_prim = new double[n];   
    double* b_prim = new double[n];   

   for (int i=0; i<n; i++)
    {
        double x = h*(i+1);
        f[i] = func_f(x);
    }


    for (int i=0; i<n-1; i++)
    {
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;
    }
    b[n-1] = 2;

    b_prim[0] = b[0];
    for (int i=0; i<n-1; i++)
    {
        b_prim[i+1] = 2 - 1/2; // (29.08.19) Are we allowed to assume b[i]=2 and precalculate b_prim[i]? (30.08.19) This is wrong!
        f_prim[i+1] = f[i+1] - 1/2;
    }
    for (int i=n-1; i>=1; i--)
    {
        f_prim[i-1] = f_prim[i-1] - 1/b_prim[i]*f_prim[i];
    }
    for (int i=0; i<n; i++)
    {
        v[i] = f_prim[i]/b_prim[i];
    }

    int value = n;
    std::string filename = "Poisson_values_n_" + std::to_string(n)
                             + ".txt";
    write_to_file(filename, n, v);

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] f;
    delete[] v;
    delete[] f_prim;
    delete[] b_prim;
}
*/

void write_to_file(std::string filename, int n, double*computed_val)
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

/*
// Old solution
void write_to_file(std::string filename, int n, double*v)
{
    double* u   = new double[n];
    double* eps = new double[n];

    double h    = 1.0/(n+1);

    for (int i=0; i<n; i++)
    {
        double x = h*(i+1);
        u[i]   = func_u(x);
        if(u[i] == 0)
        {
            eps[i] = -1;
        }
        else
        {
            eps[i] = epsilon(v[i], u[i]);
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
        results_file << std::setw(25) << std::setprecision(16) << u[i]
                     << std::setw(25) << std::setprecision(16) << v[i]
                     << std::setw(25) << std::setprecision(16) << eps[i]
                     << std::endl;
    }
    results_file << std::setw(25) << 0
                 << std::setw(25) << 0
                 << std::setw(25) << 0 << std::endl;
    results_file.close();
}
*/

double relative_error(double computed_val, double exact_val)
{
    return fabs((computed_val - exact_val)/exact_val);
}

/*
// Old solution
double epsilon(double v, double u)
{
    return (fabs((v - u)/u));
}
*/

int main(int argc, char *argv[]) 
{
    int n = atoi(argv[1]);
    gauss_elim(n);
    gauss_elim_spes(n);
 
    return 1;
}