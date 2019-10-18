#include <iostream>
#include <random>
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


int mc_integration()
{
    /*
    Monte Carlo integration of the function exp(-2*2*(r1 + r2))/|r1 - r2|.

    Returns
    -------
    seed : int
        The seed is the system time in seconds from UNIX epoch.
    */

    int N = 1e7;       // number of iterations
    
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
    double integral_sum_square = 0;
    double N5 = std::pow(N, 5);         // pre-calculated for the inner sum
    double integrand_tmp;               // temporary value for squaring the integrand so we don't have to call the function twice   
    int counter = 0;
    
    
    
    for (int i = 0; i < N; i++)
    {

        if ( counter%(N/10) == 0 )
        {
            std::cout << "iteration: " << i << " of " << N-1 << std::endl;
        }
        // drawing random numbers from the distributions
        // drawing twice for each variable for calculating the variance
        double r1 = exp_dist(engine); double R1 = exp_dist(engine);
        double r2 = exp_dist(engine); double R2 = exp_dist(engine);
        double theta1 = uniform_theta(engine); double Theta1 = uniform_theta(engine);
        double theta2 = uniform_theta(engine); double Theta2 = uniform_theta(engine);
        double phi1 = uniform_phi(engine); double Phi1 = uniform_phi(engine);
        double phi2 = uniform_phi(engine); double Phi2 = uniform_phi(engine);

        // adding to the integrand sum
        integral_sum += integrand(r1, r2, theta1, theta2, phi1, phi2)
            *r1*r1*r2*r2*std::sin(theta1)*std::sin(theta2);
        
        // adding to the f(x)**2 sum for the variance
        integrand_tmp = integrand(R1, R2, Theta1, Theta2, Phi1, Phi2);
        integral_sum_square += integrand_tmp*integrand_tmp*R1*R2*R1*R2
            *std::sin(Theta1)*std::sin(Theta2);

        counter++;
    }

    
    integral_sum /= N;                 // number of samples
    integral_sum_square /= N;          // number of samples
    
    integral_sum_square /= std::pow( (2*2), 10 );   // (2*alpha)**10
    integral_sum /= std::pow( (2*2), 5);            // (2*alpha)**5

    double variance = (integral_sum_square - integral_sum*integral_sum)*4*std::pow(pi, 4);
    
    std::cout << "\nvariance: " << variance << std::endl;
    std::cout << "std: " << std::sqrt(variance) << std::endl;

    integral_sum *= 4*std::pow(pi, 4);     // theta, phi interval
    
    std::cout << "calculated: " << integral_sum << std::endl;
    std::cout << "correct answer: " << 5*pi*pi/(16*16) << std::endl;
    std::cout << "error: " << std::fabs(integral_sum - 5*pi*pi/(16*16)) << std::endl;
    std::cout << "iterations: " << N << std::endl;

    return seed;
}


int main()
{
    mc_integration();
    return 1;
}
