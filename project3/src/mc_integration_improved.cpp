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
        The value of the integrand in the spesified point. (exp(-2*2*(r1 + r2))/|r1 - r2|)
    */
    double beta = std::cos(theta1)*std::cos(theta2) + std::sin(theta1)*std::sin(theta2)*std::cos(phi1 - phi2);
    double r12  = r1*r1 + r2*r2 - 2*r1*r2*beta;
    
    if (r12 == 0)
    {
        return 0;
    }    
    else
    {
        r12 = std::sqrt(r12);
        
        // return 1/r12;
        return std::exp(-2*2*(r1 + r2))/r12;
    }
    
}


void mc_integration_improve()
{
    /*
    Monte Carlo integration of the function exp(-2*2*(r1 + r2))/|r1 - r2|.
    */
    int N = 15; // grid points
    // integral limits, approx. infinity
    float a = 2;
    double lambda = 2; // NB! This must be spesified correctly.

    // Generate engine
    int seed = 1424;
    std::mt19937 engine(seed);

    // Generate distributions
    std::uniform_real_distribution<double> uniform(0, a);
    std::exponential_distribution<double> exp_dist(lambda);


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
                            double r0 = uniform(engine);
                            double r1 = uniform(engine);
                            double theta0 = uniform(engine);
                            double theta1 = uniform(engine);
                            double phi0 = uniform(engine);
                            double phi1 = uniform(engine);

                            // Multiplying the weights with the integrand.
                            gauss_sum += integrand(r0, r1, theta0, theta1, phi0, phi1);
                        }
                    }
                }
            }
        }
    }
    
    gauss_sum /= pow(N, 6);
    gauss_sum /= pow(1/a, 2);
    gauss_sum /= pow(1/(2*pi), 2);
    gauss_sum /= pow(1/pi, 2);

    std::cout << "calculated: " << gauss_sum << std::endl;
    std::cout << "correct answer: " << 5*pi*pi/(16*16) << std::endl;
    std::cout << "error: " << std::fabs(gauss_sum - 5*pi*pi/(16*16)) << std::endl;
}


int main()
{
    mc_integration_improve();
    return 1;
}