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

    int N = 10;       // number of iterations = N**6
    float lambda = 1; // For the exp. distribution function.

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
    double N5 = std::pow(N, 5); // pre-calculated for the inner sum

    for (int i0 = 0; i0 < N; i0++)
    {   // first loop is for displaying progress info without an if statement
        
        std::cout << "outer loop: " << i0 << " of " << N-1 << std::endl;

        for (int i1 = 0; i1 < N5; i1++)
        {
            // drawing random numbers from the distributions
            double r1 = exp_dist(engine); double R1 = exp_dist(engine);
            double r2 = exp_dist(engine); double R2 = exp_dist(engine);
            double theta1 = uniform_theta(engine); double Theta1 = uniform_theta(engine);
            double theta2 = uniform_theta(engine); double Theta2 = uniform_theta(engine);
            double phi1 = uniform_phi(engine); double Phi1 = uniform_phi(engine);
            double phi2 = uniform_phi(engine); double Phi2 = uniform_phi(engine);

            // adding to the integrand sum
            integral_sum += integrand(r1, r2, theta1, theta2, phi1, phi2)
                *r1*r1*r2*r2*std::sin(theta1)*std::sin(theta2);

            integral_sum_square += integrand(R1,R2,Theta1,Theta2,Phi1,Phi2)*integrand(R1,R2,Theta1,Theta2,Phi1,Phi2)
                *R1*R1*R2*R2*std::sin(Theta1)*std::sin(Theta2)*R1*R1*R2*R2*std::sin(Theta1)*std::sin(Theta2);        
        }
    }

    integral_sum /= std::pow(N, 6);        // number of samples
    integral_sum_square /= std::pow(N, 6);        // number of samples
    std::cout << "\nvarians: " << integral_sum_square - integral_sum*integral_sum << std::endl;
    
    integral_sum *= 4*std::pow(pi, 4);     // theta, phi interval
    integral_sum /= std::pow((2*2), 5);    // (2*alpha)**5

    std::cout << "calculated: " << integral_sum << std::endl;
    std::cout << "correct answer: " << 5*pi*pi/(16*16) << std::endl;
    std::cout << "error: " << std::fabs(integral_sum - 5*pi*pi/(16*16)) << std::endl;
    // std::cout << "iterations: " << std::pow(N, 6) << std::endl;

    return seed;
}


int main()
{
    mc_integration();
    return 1;
}
