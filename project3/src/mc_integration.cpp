#include <iostream>
#include <random>
#include <time.h>
#include <chrono>
#include <fstream>
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


void mc_integration()
{
    /*
    Monte Carlo integration of the function exp(-2*2*(r1 + r2))/|r1 - r2|.

    Returns
    -------
    seed : int
        The seed is the system time in seconds from UNIX epoch.
    */

    // generating data file
    std::ofstream mc_data_file;
    mc_data_file.open("mc_data.txt", std::ios_base::app);

    // writing title to file
    mc_data_file << std::setw(20) << "N" << std::setw(20) << "error";
    mc_data_file << std::setw(20) << "calculated";
    mc_data_file << std::setw(20) << "exact";
    mc_data_file << std::setw(20) << "comp time (s)" << std::endl;

    // int N = 10;     // number of iterations = N**6
    int N_end = 40;

    for (int N = 1; N <= N_end; N = N + 2)
    {   // loops over grid values

        // starting timer
        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
        
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

        // ending timer
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        std::chrono::duration<double> comp_time = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);

        double exact = 5*pi*pi/(16*16);
        double error = std::fabs(integral_sum - exact);

        std::cout << "\ncalculated: " << integral_sum << std::endl;
        std::cout << "correct answer: " << exact << std::endl;
        std::cout << "error: " << error << std::endl;
        std::cout << N << " of " << N_end << "\n" << std::endl;

        mc_data_file << std::setw(20) << N << std::setw(20) << error;
        mc_data_file << std::setw(20) << integral_sum;
        mc_data_file << std::setw(20) << exact;
        mc_data_file << std::setw(20) << comp_time.count() << std::endl;

    }

}


int main()
{
    mc_integration();
    return 0;
}