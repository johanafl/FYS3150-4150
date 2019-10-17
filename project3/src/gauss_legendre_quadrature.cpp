#include <cmath>
#include <iostream>
#include <chrono>
#include <fstream>

double const pi = 3.14159265359; 


void gauss_legendre_points(double x1, double x2, double x[], double w[], int n)
{
    int m, j, i;
    double z1, z, xm, xl, pp, p3, p2, p1;
    double *x_low, *x_high, *w_low, *w_high;
    double const ZERO = std::pow(10, -10);

    m  = (n + 1)/2;          // roots are symmetric in the interval
    xm = 0.5 * (x2 + x1);
    xl = 0.5 * (x2 - x1);

    x_low  = x;              // pointer initialization
    x_high = x + n - 1;
    w_low  = w;
    w_high = w + n - 1;

    for(i = 1; i <= m; i++)
    {      // loops over desired roots
        z = std::cos(pi * (i - 0.25)/(n + 0.5));

        /*
        Starting with the above approximation to the ith root
        we enter the main loop of refinement bt Newtons method.
        */

        do
        {
            p1 = 1.0;
	        p2 = 0.0;

   	        /*
	        loop up recurrence relation to get the
            Legendre polynomial evaluated at x
            */

            for(j = 1; j <= n; j++)
            {
	            p3 = p2;
	            p2 = p1;
	            p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
	        }

	        /*
	        p1 is now the desired Legrendre polynomial. Next compute
            ppp its derivative by standard relation involving also p2,
            polynomial of one lower order.
            */
 
	        pp = n*(z*p1 - p2)/(z*z - 1.0);
	        z1 = z;
	        z  = z1 - p1/pp;                   // Newton's method
        }
        while(std::fabs(z - z1) > ZERO);

        /* 
	    ** Scale the root to the desired interval and put in its symmetric
        ** counterpart. Compute the weight and its symmetric counterpart
        */

      *(x_low++)  = xm - xl*z;
      *(x_high--) = xm + xl*z;
      *w_low      = 2.0*xl/((1.0 - z*z)*pp*pp);
      *(w_high--) = *(w_low++);
    }
}


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


void gauss_legendre_quadrature()
{   /*
    Calculates an integral using Gauss-Legendre quadurature. Loops over a set of
    grid point values (N), and writes the calculated and analytical result to
    a file along with the error.
    */


    // generating data file
    std::ofstream legendre_data_file;
    legendre_data_file.open("legendre_data.txt", std::ios_base::app);

    // writing title to file
    legendre_data_file << std::setw(20) << "N" << std::setw(20) << "error";
    legendre_data_file << std::setw(20) << "calculated";
    legendre_data_file << std::setw(20) << "exact";
    legendre_data_file << std::setw(20) << "comp time (s)" << std::endl;

    // int N = 10;                 // grid points

    int N_end = 40;

    for (int N = 1; N <= N_end; N = N + 2)
    {   // loops over grid values

        // starting timer
        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
        
        double *x = new double [N]; // array of x values
        double *w = new double [N]; // array of weights

        // integral limits, approx. infinity
        float a = -2;
        float b = -a;
        
        // Finding the weights and points for integration.
        gauss_legendre_points(a, b, x, w, N);

        double integral_sum = 0;

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
                                // Multiplying the weights with the integrand.
                                integral_sum += w[i0]*w[i1]*w[i2]*w[i3]*w[i4]*w[i5]
                                    *integrand(x[i0], x[i1], x[i2], x[i3], x[i4], x[i5]);
                            }
                        }
                    }
                }
            }
        }

        // ending timer
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        std::chrono::duration<double> comp_time = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);

        double exact = 5*pi*pi/(16*16);
        double error = std::fabs(integral_sum - exact);

        std::cout << "\ncalculated: " << integral_sum << std::endl;
        std::cout << "correct answer: " << exact << std::endl;
        std::cout << "error: " << error << std::endl;
        std::cout << N << " of " << N_end << "\n" << std::endl;

        legendre_data_file << std::setw(20) << N << std::setw(20) << error;
        legendre_data_file << std::setw(20) << integral_sum;
        legendre_data_file << std::setw(20) << exact;
        legendre_data_file << std::setw(20) << comp_time.count() << std::endl;
        
        delete[] x;
        delete[] w;

    }

}




int main()
{   
    gauss_legendre_quadrature();
    return 0;
}