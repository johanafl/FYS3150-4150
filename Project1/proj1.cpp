
void set_arrays(int n) 
{
    a      = new double[n-1];
    b      = new double[n];
    c      = new double[n-1];
    f      = new double[n];
    v      = new double[n];
    f_prim = new double[n];
    b_prim = new double[n];
}

void gauss_elim(int n) 
{
    b_prim[0] = b[0];
    for (int i=1; i<=n; i++)
    {
        b_prim[i+1] = b[i+1] - a[i]/b[i]*c[i];
        f_prim[i+1] = b[i+1] - a[i]/b[i]*c[i];
    }
    for (int i=n; i>=1; i--)
    {
        f_prim[i-1] = f_prim[i-1] - c[i-1]/b_prim[i]*f_prim[i];
    }
    for (int i=0; i<=n; i++)
    {
        v[i] = f_prim[i]/b_prim[i];
    }
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
    return 0;
}