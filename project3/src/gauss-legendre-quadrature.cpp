#include <cmath>
#include <iostream>
#define pi 3.14159265359


void gauss_legendre_points(double x1, double x2, double x[], double w[], int n)
{
    int m, j, i;
    double z1, z, xm, xl, pp, p3, p2, p1;
    //double const pi = 3.14159265359; 
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
{   /*
    The integrand to be integrated.

    Parameters
    ----------
    many parameters, will write this when things seem to work
    */

    // dist is the distance between r1 and r2 vectors
    double dist = (x0 - x3)*(x0 - x3) + (x1 - x4)*(x1 - x4) + (x2 - x5)*(x2 - x5);
    
    if (dist == 0)
    {   // sets integrand to 0 if x0 = x3, x1 = x4, and x2 = x5
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
{   
    int N = 35;                 // grid points
    double *x = new double [N]; // array of x values
    double *w = new double [N]; // array of weights

    // integral limits, approx. infinity
    float a = -2;
    float b = -a;
    
    gauss_legendre_points(a, b, x, w, N);


    double gauss_sum = 0;

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
                            gauss_sum += w[i0]*w[i1]*w[i2]*w[i3]*w[i4]*w[i5]
                                *integrand(x[i0], x[i1], x[i2], x[i3], x[i4], x[i5]);
                        }
                    }
                }
            }
        }
    }

    std::cout << "calculated: " << gauss_sum << std::endl;
    std::cout << "correct answer: " << 5*pi*pi/(16*16) << std::endl;
    std::cout << "error: " << std::fabs(gauss_sum - 5*pi*pi/(16*16)) << std::endl;
    delete[] x;
    delete[] w;
}




int main()
{   
    gauss_legendre_quadrature();
    



    return 0;
}