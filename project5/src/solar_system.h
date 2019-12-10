#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include "solver.h"

const double AU = 149597871e3;  // Meters in one AU.
const double yr = 31556926;     // Seconds in a year.
const double pi = 3.14159265358979323846;
const double c  = 63197.790926112524;   // Speed of light in vacuum, [AU/yr].
// const double G_circular  = 4*pi*pi;  // Gravitational constant, [AU^3/yr^2].
const double solar_mass = 1.988e30;     // Mass of the sun, [kg].
const double G = 6.67e-11;              // Gravitational constant, [m^3/(s^2 * kg)].
const double GM = G/(AU*AU*AU)*yr*yr*solar_mass;    // [AU^3/yr^2]

class SolarSystem
/*
Keep track of all celestial objects in the solar system. Solve the solar
system.
*/
{
protected:    
    int num_planets = 0;
    int array_size  = 3;

    double* pos  = new double[array_size];
    double* vel  = new double[array_size];
    double* mass = new double[array_size/3];

    double r, x, y, z;

    // For relativistic correction to Mercury.
    arma::vec l_vec;
    double l_square, acc1, vx, vy, vz;

    // For acceleration_3.
    double rpow;
    double beta = 0;
    bool beta_set = false;

    void resize()
    {   /*
        Resize the pos, vel, mass arrays when the number of celestial
        objects increase.
        */
        array_size *= 2;
        double* tmp_pos  = new double[array_size];
        double* tmp_vel  = new double[array_size];
        double* tmp_mass = new double[array_size/3];

        for (int i = 0; i < 3*num_planets; i++)
        {   // Move position and velocity data to larger arrays.
            tmp_pos[i] = pos[i];
            tmp_vel[i] = vel[i];
        }

        for (int i = 0; i < num_planets; i++)
        {   // Move mass to larger array.
            tmp_mass[i] = mass[i];
        }

        delete[] pos;
        delete[] vel;
        delete[] mass;

        pos  = tmp_pos;
        vel  = tmp_vel;
        mass = tmp_mass;
    }

    arma::vec get_U0()
    {   /*
        Add all the initial positions and initial velocities to a
        single initial condition vector.

        Returns
        -------
        U0 : arma::vec
            Initial condition vector. r1, v1, r2, v2, ..., rn, vn.
            r = x, y, z; v = vx, vy, vz.
        */
        arma::vec U0(num_planets*3*2);
        U0.zeros();

        for (int i = 0; i < num_planets*3; i++)
        {   // Adding positions to initial condition vector.
            U0(i) = pos[i];
        }
        
        for (int i = 0; i < num_planets*3; i++)
        {   // Adding velocities to initial condition vector.
            U0(i + num_planets*3) = vel[i];
        }
        
        return U0;
    }

    arma::vec acceleration_1(arma::vec u)
    {   /*
        Acceleration for two-body problem deduced from the gravitational
        force, assuming that one object is at center at all time. That
        is, position = 0 and velocity = 0.

        Parameters
        ----------
        u : arma::vec
            Vector containing the position of the moving object in the
            x-, y- and z-direction (in that order), in astronomical
            units, [AU].

        Returns
        -------
        acc : arma::vec
            Vector containing the acceleration of the moving object in
            the x-, y- and z-direction (in that order), in astronomical
            units per years squared, [AU/yr^2].

        Note
        ----
        The object at rest is assumed to be the sun, and everything is
        therefore scaled after one solar mass.
        */

        arma::vec acc(3);
        acc.zeros();

        x = u(0);
        y = u(1);
        z = u(2);

        // Radial distance from the sun.
        r = std::sqrt(x*x + y*y + z*z);

        // Acceleration in x-, y- and z-direction.
        acc(0) -= 4*pi*pi*x/(r*r*r);
        acc(1) -= 4*pi*pi*y/(r*r*r);
        acc(2) -= 4*pi*pi*z/(r*r*r);

        return acc;
    }

    arma::vec acceleration_2(arma::vec u)
    {   /*
        Acceleration for N-body problem deduced from the gravitational
        force, assuming that one object is at center at all time. That
        is, position = 0 and velocity = 0.

        Parameters
        ----------
        u : arma::vec
            Vector containing the position of the moving objects in the
            x-, y- and z-direction, in astronomical units, [AU]. The
            vector u shoud be given as u = [x1, y1, z1, x2, y2, z2,
            ... , xN, yN, zN], where the indicies corresponds to the
            i-th object. It is important that the order is the same as
            when the object was added to the class to match up with the
            correct mass.

        Returns
        -------
        acc : arma::vec
            Vector containing the acceleration of the moving objects in
            the x-, y- and z-direction, in astronomical units per year
            squared, [AU/yr^2]. The vector acc looks like 
            acc = [ax1, ay1, az1, ax2, ay2, az2, ... , axN, ayN, azN].

        Note
        ----
        The object at rest is assumed to be the sun, and everything is
        therefore scaled after one solar mass.
        */
        arma::vec acc(3*num_planets);
        acc.zeros();

        // Acceleration due to gravitational pull from the sun.
        for (int i = 0; i < num_planets; i++)
        {
            x = u(3*i + 0);
            y = u(3*i + 1);
            z = u(3*i + 2);
            // Radial distance from the sun for the i-th object.
            r = std::sqrt(x*x + y*y + z*z);

            acc(3*i + 0) -= GM*x/(r*r*r);
            acc(3*i + 1) -= GM*y/(r*r*r);
            acc(3*i + 2) -= GM*z/(r*r*r);
        }

        // Acceleration due to gravitational pull from the j-th object
        for (int i = 0; i < num_planets; i++)
        {
            for (int j = 0; j < num_planets; j++)
            {   
                if (i == j) continue;
                
                x = u(3*j + 0) - u(3*i + 0);
                y = u(3*j + 1) - u(3*i + 1);
                z = u(3*j + 2) - u(3*i + 2);
                r = std::sqrt(x*x + y*y + z*z);

                // Acceleration in x-, y- and z-direction
                acc(3*i + 0) += GM*mass[j]*x/(r*r*r);
                acc(3*i + 1) += GM*mass[j]*y/(r*r*r);
                acc(3*i + 2) += GM*mass[j]*z/(r*r*r);
            }
        }

        return acc;
    }

    arma::vec acceleration_3(arma::vec u)
    {   /*
        The same ass acceleration_2, but with adjustable exponent of
        r (for task 5d).
        */
        
        arma::vec acc(3*num_planets);
        acc.zeros();
        // Acceleration due to gravitational pull from the sun
        for (int i = 0; i < num_planets; i++)
        {
            x = u(3*i + 0);
            y = u(3*i + 1);
            z = u(3*i + 2);
            // Radial distance from the sun for the i-th object.

            r = std::sqrt(x*x + y*y + z*z);
            rpow = std::pow(r, beta);
            
            
            acc(3*i + 0) -= GM*x/(r*rpow);
            acc(3*i + 1) -= GM*y/(r*rpow);
            acc(3*i + 2) -= GM*z/(r*rpow);
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
                    rpow = std::pow(r, beta);

                    // Acceleration in x-, y- and z-direction
                    acc(3*i + 0) += GM*mass[j]*x/(r*rpow);
                    acc(3*i + 1) += GM*mass[j]*y/(r*rpow);
                    acc(3*i + 2) += GM*mass[j]*z/(r*rpow);
                }
            }
        }
        return acc;
    }

    arma::vec acceleration_4(arma::vec u)
    {   /*
        Acceleration for N-body problem deduced from the gravitational
        force, assuming that all objects are allowed to move and the
        center of mass is at rest for all time t.

        Parameters
        ----------
        u : arma::vec
            Vector containing the position of the objects in the x-, y-
            and z-direction, in astronomical units, [AU]. The vector u
            shoud be given as
            u = [x1, y1, z1, x2, y2, z2, ... , xN, yN, zN], where the
            indicies corresponds to the i-th object. It is important
            that the order is the same as when the object was added to
            the class to match up with the correct mass.

        Returns
        -------
        acc : arma::vec
            Vector containing the acceleration of the moving objects in
            the x-, y- and z-direction. The vector acc looks like 
            acc = [ax1, ay1, az1, ax2, ay2, az2, ... , axN, ayN, azN].
        */
        arma::vec acc(3*num_planets);
        acc.zeros();
        

        // Acceleration due to gravitational pull from the j-th object.
        for (int i = 0; i < num_planets; i++)
        {
            for (int j = 0; j < num_planets; j++)
            {   
                if (i == j) continue;
                
                x = u(3*j + 0) - u(3*i + 0);
                y = u(3*j + 1) - u(3*i + 1);
                z = u(3*j + 2) - u(3*i + 2);
                r = std::sqrt(x*x + y*y + z*z);

                // Acceleration in x-, y- and z-direction.
                acc(3*i + 0) += GM*mass[j]*x/(r*r*r);
                acc(3*i + 1) += GM*mass[j]*y/(r*r*r);
                acc(3*i + 2) += GM*mass[j]*z/(r*r*r);
            }
        }

        return acc;
    }

public:
    SolarSystem() {}

    void set_beta(double beta_input)
    {   
        beta = beta_input;
        beta_set = true;
    }

    void add_celestial_body(double mass_input, arma::vec U0)
    {   /*
        Add a celestial body to the solar system. Add input initial
        conditions to arrays pos, vel and mass.

        Parameters
        ----------
        mass_input : double
            Mass of the celestial body.

        U0 : arma::vec
            Vector containing the initial position and velocity of the
            celestial body. {x, y, z, vx, vy, vz}.
        */
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

    void solve_system(int num_steps,double dt, int func_id, std::string method, 
        std::string filepath)
    {   /*
        When filename is specified, data is written to file.
        */

        solve_system(num_steps, dt, func_id, method, filepath, true);
    }

    void solve_system(int num_steps, double dt, int func_id, std::string method)
    {   /*
        When no filename is specified, no data is written to file.
        */

        std::string filepath = "unused";
        solve_system(num_steps, dt, func_id, method, filepath, false);
    }

    int solve_system(int num_steps, double dt, int func_id, std::string method,
        std::string filepath, bool write)
    {   /*
        Solve the solar system.

        Parameters
        ----------
        num_steps : int
            The number of time steps in the integration.

        dt : double
            Time step length.

        filepath : std::string
            Path to the file where data will be written.

        method : std::string
            Which integration method to use. Allowed values are
            'Forward Euler' and 'Velocity Verlet'.

        write : bool
            For toggling write to file on/off.

        func_id : int
            Choose which function to integrate.
        */

        if ((func_id == 3) and (not beta_set))
        {
            std::cout << "Beta not set! Please use set_beta(value). Exiting." << std::endl;
            return 1;
        }

        if ((func_id != 3) and (beta_set))
        {
            std::cout << "Beta is set to: " << beta << ", but the acceleration";
            std::cout << " dependent on beta is not selected. Please use func_id = 3. Exiting.";
            std::cout << std::endl;

            return 1;
        }

        if (method == "Velocity Verlet")
        {
            VelocityVerlet<SolarSystem> solved(num_steps, num_planets);
            arma::vec U0 = get_U0();
            solved.set_initial_conditions(U0);
            
            auto solve_time_1 = std::chrono::steady_clock::now();
            
            solved.solve(*this, dt, func_id);
            
            auto solve_time_2 = std::chrono::steady_clock::now();
            auto solve_time = std::chrono::duration_cast<std::chrono::duration<double> >(solve_time_2 - solve_time_1);
            std::cout << "solve time: " << solve_time.count() << " s" << std::endl;
            
            if (write)
            {
                auto write_time_1 = std::chrono::steady_clock::now();
                
                solved.write_to_file(filepath);   
                
                auto write_time_2 = std::chrono::steady_clock::now();
                auto write_time = std::chrono::duration_cast<std::chrono::duration<double> >(write_time_2 - write_time_1);
                std::cout << "write time: " << write_time.count() << " s" << std::endl;
            }
        }
        else if (method == "Forward Euler")
        {   
            ForwardEuler<SolarSystem> solved(num_steps, num_planets);
            arma::vec U0 = get_U0();
            solved.set_initial_conditions(U0);
            
            auto solve_time_1 = std::chrono::steady_clock::now();
            
            solved.solve(*this, dt, func_id);
            
            auto solve_time_2 = std::chrono::steady_clock::now();
            auto solve_time = std::chrono::duration_cast<std::chrono::duration<double> >(solve_time_2 - solve_time_1);
            std::cout << "solve time: " << solve_time.count() << " s" << std::endl;
            
            if (write)
            {
                auto write_time_1 = std::chrono::steady_clock::now();
                
                solved.write_to_file(filepath);   
                
                auto write_time_2 = std::chrono::steady_clock::now();
                auto write_time = std::chrono::duration_cast<std::chrono::duration<double> >(write_time_2 - write_time_1);
                std::cout << "write time: " << write_time.count() << " s" << std::endl;
            }
        }

        return 0;
    }

    void solve_system_mercury(int num_steps, double dt, std::string filepath)
    {   /*
        Solve the solar system.

        Parameters
        ----------
        num_steps : int
            The number of time steps in the integration.

        dt : double
            Time step length.

        filepath : std::string
            Path to the file where data will be written.

        method : std::string
            Which integration method to use. Allowed values are
            'Forward Euler' and 'Velocity Verlet'.

        write : bool
            For toggling write to file on/off.
        */

        VelocityVerlet<SolarSystem> solved(num_steps, num_planets);
        arma::vec U0 = get_U0();
        solved.set_initial_conditions(U0);
                
        solved.solve_mercury(*this, dt);
                
        solved.write_to_file(filepath);   
    }

    arma::vec acceleration(arma::vec u, double t, int func_id)
    {   /*
        The current acceleration of the system.

        Parameters
        ----------
        func_id : int
            Choose which of the accelerations to use.
        */
        
        switch (func_id)
        {
            case 1 :
                return acceleration_1(u);   // Two body.
            case 2 :
                return acceleration_2(u);   // N body, Sun at center.
            case 3 :
                return acceleration_3(u);   // N body, Adjustable beta.
            case 4 :
                return acceleration_4(u);   // N body, CM at center.
        }
        return acceleration_2(u); // Mostly to make the compiler shut up.
    }

    arma::vec acceleration_mercury(arma::vec u_pos, arma::vec u_vel, double t)
    {   /*
        Acceleration for two-body problem. Used to calculate perihelion
        of mercury i think.

        Parameters
        ----------
        u_pos : arma::vec
            Vector containing the position of the moving object in the
            x-, y- and z-direction (in that order), in astronomical
            units, [AU].

        u_vel : arma::vec
            Vector containing the velocity of the moving object in the
            x-, y- and z-direction (in that order), in astronomical
            units, [AU].

        Returns
        -------
        acc : arma::vec
            Vector containing the acceleration of the moving object in
            the x-, y- and z-direction (in that order), in astronomical
            units per years squared, [AU/yr^2].
        */
        
        arma::vec acc(3);
        acc.zeros();
        
        x = u_pos(0); vx = u_vel(0);
        y = u_pos(1); vy = u_vel(1);
        z = u_pos(2); vz = u_vel(2);

        r = std::sqrt(x*x + y*y + z*z); // Radial distance from the sun.
        l_vec = arma::cross(u_pos, u_vel);
        l_square = arma::dot(l_vec, l_vec);

        acc1  -= GM/(r*r*r)*(1 + 3*l_square/(r*r*c*c)); 
        acc(0) = acc1*x;
        acc(1) = acc1*y;
        acc(2) = acc1*z;

        return acc;
    }

    ~SolarSystem()
    {
        delete[] pos;
        delete[] vel;
        delete[] mass;
    }
};

#endif