#include <cmath>
/*
void set_arrays(int n) 
{
    a      = new double[n-1]; // = -1
    b      = new double[n];   // = 2
    c      = new double[n-1]; // = -1
    f      = new double[n];   // = func_val*h**2, where 1/(n+1)=h
    v      = new double[n];   // = What we want to find
    f_prim = new double[n];   // = func_val for calculation (values are set in fuction below)
    b_prim = new double[n];   // = b array for calculation (values are set in fuction below)

    for (int i=0; i<n; i++)
    {
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;
    }
    b[n] = 2;
}
*/
double func_f(double x)
{
    return 100*exp(-10*x);
}

void gauss_elim(int n) 
{
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
}

void gauss_elim_spes(int n) 
{
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
}

double fuck_u(double x)
{
    return 1 - (1 - exp(-10))*x - exp(-10*x);
}

double eps(double v, double u)
{
    return log10(abs((v - u)/u));
}

/*
~ArrayList() 
{
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] f;
    delete[] v;
    delete[] f_prim;
    delete[] b_prim;
}
*/

int main(int argc, char *argv[]) 
{
    n = argv[1];
    double h = 1/(n+1);

    a      = new double[n-1]; // = -1
    b      = new double[n];   // = 2
    c      = new double[n-1]; // = -1
    f      = new double[n];   // = func_val*h**2, where 1/(n+1)=h
    v      = new double[n];   // = What we want to find
    f_prim = new double[n];   // = func_val for calculation (values are set in fuction below)
    b_prim = new double[n];   // = b array for calculation (values are set in fuction below)

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

    /* Write to file before delete!!!!!!!! */
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] f;
    delete[] v;
    delete[] f_prim;
    delete[] b_prim;

    return 1;
}