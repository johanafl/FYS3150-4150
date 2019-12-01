#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include "solver.h"

const double grav_const_G = 4*pi*pi;    // [AU^3/(yr^2 * M_sun)]
const double solar_mass   = 1.988e30;   // [kg]
const double sped_of_lit  = 63197.790926112524; // [AU/yr]

class Solarsystem
{
private:
    
    int num_planets = 0;
    int array_size  = 3;
    double* pos  = new double[array_size];
    double* vel  = new double[array_size];
    double* mass = new double[array_size/3]; // [M_sun]

public:
    Solarsystem(){}    
    void add_planet(double mass_input, arma::vec U0)
    {   
        while (array_size <= 3*num_planets + 3)
        {   // Resize all arrays if they cant fit 3 (1) more values.
            resize();
        }

        pos[3*num_planets + 0] = U0(0);
        pos[3*num_planets + 1] = U0(1);
        pos[3*num_planets + 2] = U0(2);

        vel[3*num_planets + 0] = U0(3);
        vel[3*num_planets + 1] = U0(4);
        vel[3*num_planets + 2] = U0(5);

        mass[num_planets] = mass_input/solar_mass;

        num_planets++;
    }

    arma::vec get_U0()
    {   
        arma::vec U0(num_planets*3*2);
        U0.zeros();

        for (int i = 0; i < num_planets*3; i++)
        {
            U0(i) = pos[i];
        }
        
        for (int i = 0; i < num_planets*3; i++)
        {
            U0(i + num_planets*3) = vel[i];
        }
        return U0;
    }

    void resize()
    {   
        array_size *= 2;
        double* tmp_pos  = new double[array_size];
        double* tmp_vel  = new double[array_size];
        double* tmp_mass = new double[array_size/3];

        for (int i = 0; i < 3*num_planets; i++)
        {
            tmp_pos[i] = pos[i];
            tmp_vel[i] = vel[i];
        }

        for (int i = 0; i < num_planets; i++)
        {
            tmp_mass[i] = mass[i];
        }

        delete[] pos;
        delete[] vel;
        delete[] mass;

        pos  = tmp_pos;
        vel  = tmp_vel;
        mass = tmp_mass;
    }

    arma::vec acceleration_1(arma::vec u)
    {   /*
        Acceleration for two-body problem deduced from the gravitational force, 
        assuming that one object is at center at all time. That is, 
        position = 0 and velocity = 0.

        Parameters
        ----------
        u : arma::vec
            Vector containing the position of the moving object in the x-, y- 
            and z-direction (in that order), in astronomical units, [AU].

        Returns
        -------
        value : arma::vec
            Vector containing the acceleration of the moving object in the x-, 
            y- and z-direction (in that order), in astronomical units per years
            squared, [AU/yr^2].

        Note
        ----
        The object at rest is assumed to be the sun, and everything is therefore
        scaled after one solar mass.
        */
        arma::vec value(3);
        value.zeros();
        double r;
        double x;
        double y;
        double z;

        x = u(0);
        y = u(1);
        z = u(2);
        r = std::sqrt(x*x + y*y + z*z); // Radial distance from the sun

        value(0) -= 4*pi*pi*x/(r*r*r); // Acceleration in x-, y- and z-direction
        value(1) -= 4*pi*pi*y/(r*r*r);
        value(2) -= 4*pi*pi*z/(r*r*r);

        return value;
    }

    arma::vec acceleration_2(arma::vec u)
    {   /*
        Acceleration for N-body problem deduced from the gravitational force, 
        assuming that one object is at center at all time. That is, 
        position = 0 and velocity = 0.

        Parameters
        ----------
        u : arma::vec
            Vector containing the position of the moving objects in the x-, y- 
            and z-direction, in astronomical units, [AU]. The vector u shoud be 
            given as u = [x1, y1, z1, x2, y2, z2, ... , xN, yN, zN], where the
            indicies corresponds to the i-th object. It is important that the
            order is the same as when the object was added to the class to match
            up with the correct mass.

        Returns
        -------
        accel : arma::vec
            Vector containing the acceleration of the moving objects in the x-, 
            y- and z-direction, in astronomical units per years
            squared, [AU/yr^2]. The vector accel looks like 
            accel = [ax1, ay1, az1, ax2, ay2, az2, ... , axN, ayN, azN].

        Note
        ----
        The object at rest is assumed to be the sun, and everything is therefore
        scaled after one solar mass.
        */
        arma::vec accel(3*num_planets);
        accel.zeros();

        double r;
        double x;
        double y;
        double z;

        // Acceleration due to gravitational pull from the sun
        for (int i = 0; i < num_planets; i++)
        {
            x = u(3*i + 0);
            y = u(3*i + 1);
            z = u(3*i + 2);
            r = std::sqrt(x*x + y*y + z*z); // Radial distance from the sun for the i-th object.

            // unsure of how to get the sign correct for x, y and z.
            accel(3*i + 0) -= grav_const_G*x/(r*r*r);
            accel(3*i + 1) -= grav_const_G*y/(r*r*r);
            accel(3*i + 2) -= grav_const_G*z/(r*r*r);
        }

        // Acceleration due to gravitational pull from the j-th object
        for (int i = 0; i < num_planets; i++)
        {
            for (int j = 0; j < num_planets; j++)
            {
                if (i != j)
                {
                    x = u(3*j + 0) - u(3*i + 0);
                    y = u(3*j + 1) - u(3*i + 1);
                    z = u(3*j + 2) - u(3*i + 2);
                    r = std::sqrt(x*x + y*y + z*z);

                    // unsure of how to get the sign correct for x, y and z.
                    accel(3*i + 0) -= grav_const_G*mass[j]*x/(r*r*r); // Acceleration in x-, y- and z-direction
                    accel(3*i + 1) -= grav_const_G*mass[j]*y/(r*r*r);
                    accel(3*i + 2) -= grav_const_G*mass[j]*z/(r*r*r);
                }
            }
        }

        return accel;
    }

    arma::vec acceleration_3(arma::vec u)
    {   /*
        Acceleration for N-body problem deduced from the gravitational force, 
        assuming that all objects are allowed to move and the center of mass is
        at rest for all time t.

        Parameters
        ----------
        u : arma::vec
            Vector containing the position of the objects in the x-, y- and 
            z-direction, in astronomical units, [AU]. The vector u shoud be 
            given as u = [x1, y1, z1, x2, y2, z2, ... , xN, yN, zN], where the
            indicies corresponds to the i-th object. It is important that the
            order is the same as when the object was added to the class to match
            up with the correct mass.

        Returns
        -------
        accel : arma::vec
            Vector containing the acceleration of the moving objects in the x-, 
            y- and z-direction, in astronomical units per years
            squared, [AU/yr^2]. The vector accel looks like 
            accel = [ax1, ay1, az1, ax2, ay2, az2, ... , axN, ayN, azN].

        Note
        ----
        Everything is scaled after one solar mass.
        */
        arma::vec accel(3*num_planets);
        accel.zeros();

        double r;
        double x;
        double y;
        double z;

        // Acceleration due to gravitational pull from the j-th object
        for (int i = 0; i < num_planets; i++)
        {
            for (int j = 0; j < num_planets; j++)
            {
                x = u(3*j + 0) - u(3*i + 0);
                y = u(3*j + 1) - u(3*i + 1);
                z = u(3*j + 2) - u(3*i + 2);
                r = std::sqrt(x*x + y*y + z*z); // Radial distance from the sun for the i-th object.

                // unsure of how to get the sign correct for x, y and z.
                accel(3*i + 0) -= grav_const_G*mass[j]*x/(r*r*r); // Acceleration in x-, y- and z-direction
                accel(3*i + 1) -= grav_const_G*mass[j]*y/(r*r*r);
                accel(3*i + 2) -= grav_const_G*mass[j]*z/(r*r*r);
            }
        }

        return accel;
    }

    arma::vec acc_mercury(arma::vec u_pos, arma::vec u_vel)
    {   /*
        Acceleration for two-body problem. Used to calculate perihelion of 
        mercury i think.

        Parameters
        ----------
        u_pos : arma::vec
            Vector containing the position of the moving object in the x-, y- 
            and z-direction (in that order), in astronomical units, [AU].

        u_vel : arma::vec
            Vector containing the velocity of the moving object in the x-, y- 
            and z-direction (in that order), in astronomical units, [AU].

        Returns
        -------
        value : arma::vec
            Vector containing the acceleration of the moving object in the x-, 
            y- and z-direction (in that order), in astronomical units per years
            squared, [AU/yr^2].
        */
        arma::vec value(3);
        value.zeros();
        double r; arma::vec l_vec; double l_square;
        double value1;
        double x; double vx;
        double y; double vy;
        double z; double vz;

        x = u_pos(0); vx = u_vel(0);
        y = u_pos(1); vy = u_vel(1);
        z = u_pos(2); vz = u_vel(2);
        r = std::sqrt(x*x + y*y + z*z); // Radial distance from the sun
        l_vec = arma::cross(u_pos, u_vel);
        l_square = arma::dot(l_vec, l_vec);

        value1 -= 4*pi*pi/(r*r*r)*(1 + 3*l_square/(r*r*sped_of_lit*sped_of_lit)); 
        value(0) = value1*x;
        value(1) = value1*y;
        value(2) = value1*z;

        return value;
    }

    arma::vec acceleration(arma::vec u, double t)
    {
        return acceleration_2(u);
    }

    void solve_system()
    {   
        int num_steps = 1e5;
        double dt = 1e-3;
        arma::vec U0 = get_U0();
        
        // std::cout << num_planets << std::endl;  
        // ForwardEuler<Solarsystem> solved(num_steps, num_planets);
        VelocityVerlet<Solarsystem> solved(num_steps, num_planets);
        solved.set_initial_conditions(U0);
        // U0.print();
        // solved.pos.col(0).print();
        solved.solve(*this, dt);
        solved.write_to_file();   
    }

};

#endif