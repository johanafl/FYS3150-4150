#include <iostream>
#include <random>
#include <time.h>
double const pi = 3.14159265359; 



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
    {   // avoids division by zero
        return 0;
    } 
    
    else
    {
        double res; 
        res  = std::sqrt(x0*x0 + x1*x1 + x2*x2) + std::sqrt(x3*x3 + x4*x4 + x5*x5);
        res  = std::exp(-2*2*(res));
        res /= std::sqrt(dist);
        
        return res;
    }

}


int mc_integration()
{
    /*
    Monte Carlo integration of the function exp(-2*2*(r1 + r2))/|r1 - r2|.

    Returns
    -------
    seed : int
        The seed is the system time in seconds from UNIX epoch.
    */

    int N = 15;     // number of iterations = N**6
    
    // integral limits, approx. infinity
    float a = -2;
    float b = -a;

    // generate engine with pseudo-random seed taken from system time
    time_t seed;
    time(&seed);
    std::mt19937 engine(seed);

    // generating distribution
    std::uniform_real_distribution<double> uniform(a, b);


    double integral_sum = 0;
    double N5 = std::pow(N, 5);

    for (int i0 = 0; i0 < N; i0++)
    {   // first loop is for displaying progress info without an if statement
    
        std::cout << "outer loop: " << i0 << " of " << N-1 << std::endl;

        for (int i1 = 0; i1 < N5; i1++)
        {
            // drawing random numbers from the uniform distribution
            double x0 = uniform(engine);
            double x1 = uniform(engine);
            double x2 = uniform(engine);
            double x3 = uniform(engine);
            double x4 = uniform(engine);
            double x5 = uniform(engine);

            // adding to the integrand sum
            integral_sum += integrand(x0, x1, x2, x3, x4, x5);   
        
        }
    }


    integral_sum /= std::pow(N, 6);     // number of samples
    integral_sum *= pow((b - a), 6);    // integral interval

    std::cout << "calculated: " << integral_sum << std::endl;
    std::cout << "correct answer: " << 5*pi*pi/(16*16) << std::endl;
    std::cout << "error: " << std::fabs(integral_sum - 5*pi*pi/(16*16)) << std::endl;

    return seed;
}


int main()
{
    mc_integration();
    return 0;
}