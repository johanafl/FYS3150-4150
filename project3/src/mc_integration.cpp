#include <iostream>
#include <random>
double const pi = 3.14159265359; 

// // Generate engine
// int seed = 1424;
// mt19937 engine(seed);

// // Generate distributions
// uniform_real_distribution<float> rand_float(0, 2);
// normal_distribution<float> rand_normal(0, 3);


double integrand(double x0, double x1, double x2, double x3, double x4, double x5)
{   
    /*
    Here we calculate the integrand that we want to integrate.

    Parameters
    ----------
    x0 : double
        The x-value for the first part of the integrand/for the first electron.

    x1 : double
        The y-value for the first part of the integrand/for the first electron.

    x2 : double
        The z-value for the first part of the integrand/for the first electron.

    x3 : double
        The x-value for the second part of the integrand/for the second electron.

    x4 : double
        The y-value for the second part of the integrand/for the second electron.

    x5 : double
        The z-value for the second part of the integrand/for the second electron.


    Returns
    -------
    : double
        The value of the integrand in the spesified point. (exp(-2*2*(r1 + r2))/|r1 - r2|)
    */

    // dist is the distance between r1 and r2 vectors
    double dist = (x0 - x3)*(x0 - x3) + (x1 - x4)*(x1 - x4) + (x2 - x5)*(x2 - x5);
    
    if (dist == 0)
    {   
        // sets integrand to 0 if x0 = x3, x1 = x4, and x2 = x5
        return 0;
    }
    
    else
    {
        double res = std::sqrt(x0*x0 + x1*x1 + x2*x2) + std::sqrt(x3*x3 + x4*x4 + x5*x5);
        res = std::exp(-2*2*(res));
        res /= std::sqrt(dist);
        
        return res;
    }

}


void mc_integration()
{
    /*
    Monte Carlo integration of the function exp(-2*2*(r1 + r2))/|r1 - r2|.
    */
    int N = 15; // grid points
    // integral limits, approx. infinity
    float a = -2;
    float b = -a;

    // Generate engine
    int seed = 1424;
    std::mt19937 engine(seed);

    // Generate distributions
    std::uniform_real_distribution<double> uniform(a, b);


    double gauss_sum = 0;

    // The actual integral is approximated with a sum
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
                            double x0 = uniform(engine);
                            double x1 = uniform(engine);
                            double x2 = uniform(engine);
                            double x3 = uniform(engine);
                            double x4 = uniform(engine);
                            double x5 = uniform(engine);

                            // Multiplying the weights with the integrand.
                            gauss_sum += integrand(x0, x1, x2, x3, x4, x5);
                        }
                    }
                }
            }
        }
    }
    
    gauss_sum /= pow(N, 6);
    gauss_sum /= pow(1/(b - a), 6);

    std::cout << "calculated: " << gauss_sum << std::endl;
    std::cout << "correct answer: " << 5*pi*pi/(16*16) << std::endl;
    std::cout << "error: " << std::fabs(gauss_sum - 5*pi*pi/(16*16)) << std::endl;
}


int main()
{
    mc_integration();
    return 1;
}