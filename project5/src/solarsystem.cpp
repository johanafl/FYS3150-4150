#include "solver.h"
#include "planet.h"
#include <armadillo>

class Solarsystem
{
private:
    
    int num_planets = 0;
    int array_size  = 2;
    Planet* planets = new Planet[array_size];

    // double* ax = new double[array_size];
    // double* ay = new double[array_size];
    // double* az = new double[array_size];

    arma::vec value = arma::zeros<arma::vec>(array_size);

    Solarsystem() {}

    void add_planet(double mass, double pos_x, double pos_y, double pos_z, double vel_x, double vel_y, double vel_z)
    {
        Planet new_planet(mass, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z);
        
        // resize
        if (num_planets == array_size)
        {
            array_size *= 2;
            Planet* tmp = new Planet[array_size];

            // double* tmp_ax = new double[array_size];
            // double* tmp_ay = new double[array_size];
            // double* tmp_az = new double[array_size];

            for (int i=0; i<num_planets; i++)
            {
                tmp[i] = planets[i];
            }

            delete[] planets;
            // delete[] ax;
            // delete[] ay;
            // delete[] az;

            planets = tmp;

            // ax = tmp_ax;
            // ay = tmp_ay;
            // az = tmp_az;
        }

        planets[num_planets] = new_planet;
        num_planets += 1;
        value = arma::zeros<arma::vec>(num_planets);
    }

    arma::vec acceleration(arma::vec u, double t)
    // acceleration(double u, double t)
    {
        double G = 1;
        double r;
        double x;
        double y;
        double z;

        // arma::vec<double> value[];

        // for planet in planets:
        for (int i=0; i<num_planets; i++)
        {
            for (int j=0; j<num_planets; j++)
            {
                x = u(3*j) - u(3*i);
                y = u(3*j + 1) - u(3*i + 1);
                z = u(3*j + 2) - u(3*i + 2);
                r = std::sqrt(x*x + y*y + z*z);
                // unsure of how to get the sign correct for x, y and z.

                value(3*i)     -= G*planets[j].mass*x/(r*r*r);
                value(3*i + 1) -= G*planets[j].mass*y/(r*r*r);
                value(3*i + 2) -= G*planets[j].mass*z/(r*r*r);

                // // a += G*planet.mass/|planet.r - current_planet.r|^2
                // // Planet MUST contain position to make this work!
                // r = (planet[i].r - planet[j].r);
                // x = (planet[i].x - planet[j].x);
                // y = (planet[i].y - planet[j].y);
                // z = (planet[i].z - planet[j].z);

                // ax[j] -= G*planets[i].mass*x/(r*r*r);
                // ay[j] -= G*planets[i].mass*y/(r*r*r);
                // az[j] -= G*planets[i].mass*z/(r*r*r);
            }
        }
        return value;
    }

    
};

int main()
{
    return 0;
}