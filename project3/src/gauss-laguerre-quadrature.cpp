#include <cmath>
#include "gauss-laguerre.cpp"

double integrand(double r1, double r2, double theta1, double theta2, double phi1, double phi2)
{
    double beta = std::cos(theta1)*std::cos(theta2) + std::sin(theta1)*std::sin(theta2)*std::cos(phi1 - phi2);
    double r12  = r1*r1 + r2*r2 - 2*r1*r2*beta;
    
    if (r12 == 0)
    {
        return 0;
    }
    
    else
    {
        r12 = std::sqrt(r12);
        
        return std::exp(-2*2*(r1 + r2))/r12;
    }
    
}


void gauss_laguerre_quadrature()
{

    int N = 10;
    double alpha = 2;

    double *w[N];
    double *x[N];

    


    double gauss_sum = 0;

    for (int i0 = 0; i0 < N; i0++)
    {   
        std::cout << "outer loop: " << i0 << " of " << N-1 << std::endl;
        
        for (int i1 = 0; i1 < N; i1++)
        {
            for (int i2 = 0; i2 < N; i2++)
            {
                for (int i3 = 0; i3 < N; i3++)
                {
                    for (int i4 = 0; i4 < N; i4++)
                    {
                        for (int i5 = 0; i5 < N; i5++)
                        {
                            gauss_sum += w[i0]*w[i1]*w[i2]*w[i3]*w[i4]*w[i5]
                                *integrand(x[i0], x[i1], x[i2], x[i3], x[i4], x[i5]);
                        }
                    }
                }
            }
        }
    }

}



int main()
{   
    gauss_laguerre_quadrature();
    return 0;
}