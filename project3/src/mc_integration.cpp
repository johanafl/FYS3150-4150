#include <iostream>
#include <random>
#include <time.h>
#include <iomanip>
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

double mc_integration(int N, float lambda)
{
    /*
    Monte Carlo integration of the function exp(-2*2*(r1 + r2))/|r1 - r2|.
    */
    
    // generate engine with pseudo-random seed taken from system time
    time_t seed;
    time(&seed);
    std::mt19937 engine(seed);

    // integral limits
    float a = -lambda;
    float b = lambda;

    // generating distribution
    std::uniform_real_distribution<double> uniform(a, b);


    double integral_sum = 0;
    double integral_sum_square = 0;
    double integrand_tmp = 0;     // placeholder for temporary integrand values
    
    // double x0; double y0;
    // double x1; double y1;
    // double x2; double y2;
    // double x3; double y3;
    // double x4; double y4;
    // double x5; double y5;


    for (int i = 0; i < N+1; i++)
    {   
        if ( i%(N/10) == 0 )
        {
        std::cout << "iteration: " << i << " of " << N << std::endl;
        }
        
        // drawing random numbers from the uniform distribution
        double x0 = uniform(engine); double y0 = uniform(engine);
        double x1 = uniform(engine); double y1 = uniform(engine);
        double x2 = uniform(engine); double y2 = uniform(engine);
        double x3 = uniform(engine); double y3 = uniform(engine);
        double x4 = uniform(engine); double y4 = uniform(engine);
        double x5 = uniform(engine); double y5 = uniform(engine);

        // adding to the integrand sum
        integral_sum += integrand(x0, x1, x2, x3, x4, x5);
        // adding to the variance sum
        integrand_tmp = integrand(y0, y1, y2, y3, y4, y5);
        integral_sum_square += integrand_tmp*integrand_tmp;
    }


    integral_sum /= N;             // number of samples
    integral_sum_square /= N;      // number of samples

    double variance = (integral_sum_square - integral_sum*integral_sum)*pow( (b - a), 6);
    std::cout << "\nvariance: " << variance << std::endl;
    std::cout << "std: " << std::sqrt(variance) << std::endl;

    integral_sum *= pow( (b - a), 6);           // integral interval

    std::cout << "calculated: " << integral_sum << std::endl;
    std::cout << "correct answer: " << 5*pi*pi/(16*16) << std::endl;
    std::cout << "error: " << std::fabs(integral_sum - 5*pi*pi/(16*16)) << std::endl;

    return integral_sum;
}

class MCIntegration
{
private:
    int N_end   = 1e7;
    int dN      = 1e5;
    int N_start = dN;

    float lambda_start = 1;
    float lambda_end   = 4;
    float dlambda      = 0.2;

    bool debug = true;
    bool write_contour_data = true;
    double exact = 5*pi*pi/(16*16);

    std::ofstream mc_data_file;
    std::ofstream mc_contour_data_file;


public:
    MCIntegration()
    {
        // generating data files
        mc_data_file.open("data_files/mc_data.txt", std::ios_base::app);
        mc_contour_data_file.open("data_files/mc_contour_data.txt", std::ios_base::app);
    
    }
    
    void lambda_loop()
    {   /*
        Loops over integral limits a = -lambda, b = lambda. Approximation of
        +-infty.
        */

        mc_contour_data_file << std::setw(15) << "N range";
        mc_contour_data_file << std::setw(15) << "lambda range" << std::endl;
        mc_contour_data_file << std::setw(10) << N_start << " " << N_end << " " << dN;
        mc_contour_data_file << std::setw(10) << lambda_start << " " << lambda_end << " " << dlambda << std::endl;

        for (float lambda_current = lambda_start; lambda_current < lambda_end; lambda_current += dlambda)
        {   // loops over integral limits / infty approximations
            grid_loop(lambda_current);
        }
    }

    void grid_loop(float lambda_current, int N)
    {   /*
        Takes a single lambda and N input and computes a single run with these
        values.

        Parameters
        ----------
        lambda_current : float
            Infinity approximation value.

        N : int
            Grid point value.
        */
       
        N_start = N;
        N_end = N;
        write_contour_data = false;

        grid_loop(lambda_current);
    }

    void grid_loop(float lambda, int N_start_input, int N_end_input, int dN_input)
    {   /*
        Takes a single lambda value, and loops it over a set of grid point
        values defined by the input values.

        Parameters
        ----------
        lambda : float
            Infinity approximation value.

        N_start_input : int
            Start grid point value.

        N_end_input : int
            End grid point value.

        dN_input : int
            Grid point step size.
        */

        write_contour_data = false;

        N_start = N_start_input;
        N_end   = N_end_input;
        dN      = dN_input;
        
        grid_loop(lambda);

    }

    void grid_loop(float lambda)
    {   /*
        Loops over grid point values.

        Parameters
        ----------
        lambda : float
            Infinity approximation value.
        */


        // writing title to file
        mc_data_file << std::setw(20) << "N" << std::setw(20) << "error";
        mc_data_file << std::setw(20) << "calculated";
        mc_data_file << std::setw(20) << "exact";
        mc_data_file << std::setw(20) << "comp time (s)" << std::endl;

        double integral_sum;


        for (int N = N_start; N <= N_end; N += dN)
        {   // loops over grid values

            // starting timer
            std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
            
            // integrating
            integral_sum = mc_integration(N, lambda);
            
            // ending timer
            std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
            std::chrono::duration<double> comp_time = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
            
            
            double error = std::fabs(integral_sum - exact);

            if (debug)
            {
                // printing calculation data
                std::cout << "\ncalculated: " << integral_sum << std::endl;
                std::cout << "correct answer: " << exact << std::endl;
                std::cout << "error: " << error << std::endl;
                std::cout << "N: " <<  N << " of " << N_end << "\n";
                std::cout << "lambda: " << lambda << "\n" << std::endl;
            }


            // writing calculation data to file
            mc_data_file << std::setw(20) << N << std::setw(20) << error;
            mc_data_file << std::setw(20) << integral_sum;
            mc_data_file << std::setw(20) << exact;
            mc_data_file << std::setw(20) << comp_time.count() << std::endl;

            if (write_contour_data)
            {
                mc_contour_data_file << error << " ";
            }

        }

        if (write_contour_data)
        {
            mc_contour_data_file << std::endl;
        }

    }
};





int main()
{   
    int N = 1e6;
    float lambda = 2;

    MCIntegration q;
    q.grid_loop(lambda);

    // mc_integration(N, lambda);
    return 0;
}