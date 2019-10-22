#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include <chrono>

double const pi = 3.14159265359; 


double integrand(double r1, double r2, double theta1, double theta2, double phi1, double phi2)
{
    /*
    Here we calculate the integrand that we want to integrate.

    Parameters
    ----------
    r1 : double
        The r-value for the first part of the integrand/for the first electron.

    r2 : double
        The r-value for the second part of the integrand/for the second electron.

    theta1 : double
        The theta-value for the first part of the integrand/for the first electron.

    theta2 : double
        The theta-value for the second part of the integrand/for the second electron.

    phi1 : double
        The phi-value for the first part of the integrand/for the first electron.

    phi2 : double
        The phi-value for the second part of the integrand/for the second electron.


    Returns
    -------
    : double
        The value of the integrand in the specified point.
    */

    double tol = 1e-10;
    double cos_beta = std::cos(theta1)*std::cos(theta2) + std::sin(theta1)*std::sin(theta2)*std::cos(phi1 - phi2);
    double r12 = r1*r1 + r2*r2 - 2*r1*r2*cos_beta;
    
    if (r12 < tol)
    {   // avoids division by zero
        return 0;
    }    
    else
    {
        return 1/std::sqrt(r12);
    }
    
}


void mc_integration(int N_start, int N_end, int dN)
{
    /*
    Monte Carlo integration of the function exp(-2*2*(r1 + r2))/|r1 - r2|.

    Parameters
    ----------
    N_start : int
        Start iteration value.

    N_end : int
        End iteration value.

    dN : int
        Iteration step size.
    */

    // generating data file
    std::ofstream mc_improved_data_file;
    mc_improved_data_file.open("data_files/mc_improved_variance_data.txt", std::ios_base::app);

    // writing title to file
    mc_improved_data_file << std::setw(20) << "N" << std::setw(20) << "error";
    mc_improved_data_file << std::setw(20) << "calculated";
    mc_improved_data_file << std::setw(20) << "exact";
    mc_improved_data_file << std::setw(20) << "comp time (s)";
    mc_improved_data_file << std::setw(20) << "variance" << std::endl;

    // if we choose lambda = 1, we need only one exp distribution
    // if lambda = 2*alpha, we need another distribution with lambda = 4*alpha
    float lambda = 1;

    // generate engine with pseudo-random seed taken from system time
    time_t seed;
    time(&seed);
    std::mt19937 engine(seed);

    // generating distributions
    std::uniform_real_distribution<double> uniform_theta(0, pi);
    std::uniform_real_distribution<double> uniform_phi(0, 2*pi);
    std::exponential_distribution<double> exp_dist(lambda);

    double integral_sum = 0;
    double average_sum;
    double average_time;
    double integral_sum_square = 0;
    double integrand_tmp;               // temporary value for squaring the integrand so we don't have to call the function twice
    int average_runs = 3;      // number of iterations for the average
    double variance;
    double exact = 5*pi*pi/(16*16);

    double r1;     double R1;
    double r2;     double R2;
    double theta1; double Theta1;
    double theta2; double Theta2;
    double phi1;   double Phi1;
    double phi2;   double Phi2;

    
    for (int N = N_start; N <= N_end; N += dN)
    {   // loops over MC iterations

        // resetting the averages
        average_sum = 0;
        average_time = 0;

        for (int _ = 0; _ < average_runs; _++)
        {   // averaging to get results without too many statistical flukes

            // starting timer
            std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

            for (int i = 0; i < N; i++)
            {
                // drawing random numbers from the distributions
                // drawing twice for each variable for calculating the variance
                r1 = exp_dist(engine);          R1 = exp_dist(engine);
                r2 = exp_dist(engine);          R2 = exp_dist(engine);
                theta1 = uniform_theta(engine); Theta1 = uniform_theta(engine);
                theta2 = uniform_theta(engine); Theta2 = uniform_theta(engine);
                phi1 = uniform_phi(engine);     Phi1 = uniform_phi(engine);
                phi2 = uniform_phi(engine);     Phi2 = uniform_phi(engine);

                // adding to the integrand sum
                integral_sum += integrand(r1, r2, theta1, theta2, phi1, phi2)
                    *r1*r1*r2*r2*std::sin(theta1)*std::sin(theta2);
                
                // adding to the f(x)**2 sum for the variance
                integrand_tmp = integrand(R1, R2, Theta1, Theta2, Phi1, Phi2);
                integral_sum_square += integrand_tmp*integrand_tmp*R1*R2*R1*R2
                    *std::sin(Theta1)*std::sin(Theta2);
            }

            integral_sum /= N;                  // number of samples
            integral_sum_square /= N;          // number of samples
            integral_sum_square /= std::pow( (2*2), 10 );   // (2*alpha)**10
            integral_sum /= std::pow( (2*2), 5);            // (2*alpha)**5
            variance = (integral_sum_square - integral_sum*integral_sum)*4*std::pow(pi, 4)*4*std::pow(pi, 4);
            integral_sum *= 4*std::pow(pi, 4);     // theta, phi interval

            average_sum += integral_sum;

            // ending timer
            std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
            std::chrono::duration<double> comp_time  = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
            average_time += comp_time.count();
            
        }
        average_sum  /= average_runs;
        average_time /= average_runs;

        double error = std::fabs(average_sum - exact);

        std::cout << "\nvariance: " << variance << std::endl;
        std::cout << "std: " << std::sqrt(variance) << std::endl;
        std::cout << "calculated: " << average_sum << std::endl;
        std::cout << "correct answer: " << exact << std::endl;
        std::cout << "error: " << error << std::endl;
        std::cout << "iterations: " << N << std::endl;

        // writing calculation data to file
        mc_improved_data_file << std::setw(20) << N << std::setw(20) << error;
        mc_improved_data_file << std::setw(20) << average_sum;
        mc_improved_data_file << std::setw(20) << exact;
        mc_improved_data_file << std::setw(20) << average_time;
        mc_improved_data_file << std::setw(20) << variance << std::endl;
    }

    mc_improved_data_file.close();
}


int main()
{

    int N_end   = 1e7;
    int dN      = 5e5;
    int N_start = dN;
    
    mc_integration(N_start, N_end, dN);
    
    return 0;
}
