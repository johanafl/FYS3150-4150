
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
    //
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