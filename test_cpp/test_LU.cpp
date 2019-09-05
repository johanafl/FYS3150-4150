#include <armadillo>
#include <iostream>

int main()
{
    int n_rows = 10;
    int n_cols = 10;
    arma::mat A(n_rows, n_cols);
    for(int i=0; i<n_rows; i++)
    {
        for(int j=0; j<n_cols; j++){
            if(i == j)
            {
                A(i, j) = 2.0;
            }
            else if(i == j + 1)
            {
                A(i, j) = 3.0;
            }
            else if(i == j - 1)
            {
                A(i, j) = 1.0;
            }
        }
    }
    A.print();
    // std::cout << A.print() << std::endl ;
    return 1;
}