#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>

double func_f(double x)
{
    return 100*exp(-10*x);
}

double func_u(double x)
{
    return 1 - (1 - exp(-10))*x - exp(-10*x);
}

void write_to_file(std::string filename, int n, double*v);

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
        f_prim[i+1] = f[i+1] - a[i]/b[i]*c[i];
    }
    for (int i=n-1; i>=1; i--)
    {
        f_prim[i-1] = f_prim[i-1] - c[i-1]/b_prim[i]*f_prim[i];
    }
    for (int i=0; i<n; i++)
    {
        v[i] = f_prim[i]/b_prim[i];
    }
    
    std::string filename = "Johan_er_kul_n=" + std::to_string(n)
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
        b_prim[i+1] = 2 - 1/2; // Are we allowed to assume b[i]=2 and precalculate b_prim[i]?
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
    std::string filename = "Johan_er_kul_n=" + std::to_string(n)
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

double epsilon(double v, double u);


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

int main(int argc, char *argv[]) 
{
    int n = atoi(argv[1]);
    gauss_elim(n);
    return 1;
}

double epsilon(double v, double u)
{
    return (abs((v - u)/u));
}