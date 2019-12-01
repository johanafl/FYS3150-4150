#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include "solver.h"

class Solarsystem
{
private:
    
    int num_planets = 0;
    int array_size  = 3;
    double* pos  = new double[array_size];
    double* vel  = new double[array_size];
    double* mass = new double[array_size/3];

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

        mass[num_planets] = mass_input;

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

    arma::vec acceleration(arma::vec u, double t)
    {
        // std::cout << "hei " << std::endl;
        // double G = 1;

        arma::vec value(3);
        value.zeros();
        double r;
        double x;
        double y;
        double z;

        x = u(0);
        y = u(1);
        z = u(2);
        r = std::sqrt(x*x + y*y + z*z);

        value(0) -= 4*pi*pi*x/(r*r*r);
        value(1) -= 4*pi*pi*y/(r*r*r);
        value(2) -= 4*pi*pi*z/(r*r*r);


        return value;
    }

    void solve_system()
    {   
        int num_steps = 1e5;
        double dt = 1e-3;
        arma::vec U0 = get_U0();
        
        // std::cout << num_planets << std::endl;  
        ForwardEuler<Solarsystem> solved(num_steps, num_planets);
        // VelocityVerlet solved(num_steps);
        solved.set_initial_conditions(U0);
        // U0.print();
        // solved.pos.col(0).print();
        solved.solve(*this, dt);
        solved.write_to_file();   
    }

};

#endif