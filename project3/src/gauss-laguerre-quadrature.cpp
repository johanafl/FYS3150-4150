#include <cmath>
#include "gauss-laguerre.cpp"

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

    double tol = 1e-10;
    double cos_beta = std::cos(theta1)*std::cos(theta2) + std::sin(theta1)*std::sin(theta2)*std::cos(phi1 - phi2);
    double r12  = r1*r1 + r2*r2 - 2*r1*r2*cos_beta;
    
    if (r12 < tol)
    {
        return 0;
    }    
    else
    {
        r12 = std::sqrt(r12);
        return std::exp(-(r1 + r2) )/r12;
    }
    
}


void gauss_laguerre_quadrature()
{
    /*
    Calculate the integral of exp(-2*alpha*(r1 + r2))/|r1 - r2|, with alpha=1,
    using gaussian quadrature with laguerre polynomials for the r-dependence
    and legendre polynomials for the theta- and phi-dependence.
    */
    
    int N = 26;         // # of integration points for a single integral.
    double alpha = 2;   // laguerre assumes a function of the form x^{alpha} exp(-x), and we must specify alpha.

    double *r     = new double[N+1];  // Arrays for the r, theta and phi points.
    double *theta = new double[N];
    double *phi   = new double[N];

    double *w_r     = new double[N+1]; // Arrays for the r, theta and phi weights.
    double *w_theta = new double[N];
    double *w_phi   = new double[N];

    // Finding the weights and points for integration.
    gauss_laguerre(r, w_r, N+1, alpha);
    gauss_legendre_points(0, pi, theta, w_theta, N);
    gauss_legendre_points(0, 2*pi, phi, w_phi, N);


    double gauss_sum = 0;

    // The actual integral is approximated with a sum.
    for (int i0 = 1; i0 < N+1; i0++)
    {   
        std::cout << "outer loop: " << i0 << " of " << N << std::endl;
        
        for (int i1 = 1; i1 < N+1; i1++)
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
                            gauss_sum += w_r[i0]*w_r[i1]*w_theta[i2]*w_theta[i3]*w_phi[i4]*w_phi[i5]
                                *integrand(r[i0], r[i1], theta[i2], theta[i3], phi[i4], phi[i5])
                                *r[i0]*r[i0]*r[i1]*r[i1]*std::sin(theta[i2])*std::sin(theta[i3]);
                        }
                    }
                }
            }
        }
    }

    // factor from change of variables
    gauss_sum /= std::pow( (2*2), 5);

    std::cout << "calculated: " << gauss_sum << std::endl;
    std::cout << "correct answer: " << 5*pi*pi/(16*16) << std::endl;
    std::cout << "error: " << std::fabs(gauss_sum - 5*pi*pi/(16*16)) << std::endl;

    delete[] r;
    delete[] theta;
    delete[] phi;

    delete[] w_r;
    delete[] w_theta;
    delete[] w_phi;
}


int main()
{   
    gauss_laguerre_quadrature();
    return 1;
}