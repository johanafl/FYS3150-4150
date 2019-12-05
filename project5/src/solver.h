#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <armadillo>
#include <chrono>

template <class T>
class Solver
{
public:
    int num_steps;
    int num_stellar_objects;
    
    double dt;

    double* pos_x;
    double* pos_y;
    double* vel_x;
    double* vel_y;

    arma::mat vel;
    arma::mat pos;

    Solver(int num_steps_input, int num_stellar_objects_input) 
    : pos(3*num_stellar_objects_input, num_steps_input+1),
    vel(3*num_stellar_objects_input, num_steps_input+1)
    {
        // pos.zeros();
        // vel.zeros();
        num_steps = num_steps_input;
        num_stellar_objects = num_stellar_objects_input;
        // std::cout << num_stellar_objects << std::endl;
    }

    // Virtual to allow a method in the superclass to call a method advance
    // in a subclass.
    
    void set_initial_conditions(arma::vec U0)
    {
        // std::cout << "hei" << std::endl;
        for (int i = 0; i < 3*num_stellar_objects; i++)
        {
            pos(i, 0) = U0(i);
            // std::cout << "hei" << std::endl;
        }

        for (int i = 0; i < 3*num_stellar_objects; i++)
        {
            vel(i, 0) = U0(i + 3*num_stellar_objects);
        }
    }
    
    void solve(T object, double dt_input)
    {
        dt = dt_input;
        double t = 0;  // Currently unused.
        
        for (int k = 0; k < num_steps; k++)
        {
            advance(object, k);
        }
    }
    
    virtual void advance(T object, int k)
    {
        // Dummy method.
        std::cout << "NotImplementedError" << std::endl;
    }

    void write_to_file()
    {
        // defining data files
        std::ofstream tull_mc_tull;
        // tull_mc_tull.open("data_files/tull_mc_tull.txt", std::ios_base::app);
        tull_mc_tull.open("data_files/tull_mc_tull.txt");

        for (int i = 0; i < num_steps; i++)
        {
            tull_mc_tull << std::setw(20) << std::setprecision(15);
            tull_mc_tull << dt*i;

            for (int j = 0; j < num_stellar_objects; j++)
            {
                tull_mc_tull << std::setw(30) << std::setprecision(15);
                tull_mc_tull << pos(3*j + 0, i);
                tull_mc_tull << std::setw(30) << std::setprecision(15);
                tull_mc_tull << pos(3*j + 1, i);
                tull_mc_tull << std::setw(30) << std::setprecision(15);
                tull_mc_tull << pos(3*j + 2, i);

                tull_mc_tull << std::setw(30) << std::setprecision(15);
                tull_mc_tull << vel(3*j + 0, i);
                tull_mc_tull << std::setw(30) << std::setprecision(15);
                tull_mc_tull << vel(3*j + 1, i);
                tull_mc_tull << std::setw(30) << std::setprecision(15);
                tull_mc_tull << vel(3*j + 2, i);
            }
            tull_mc_tull << std::endl;
        }
        tull_mc_tull.close();
    }
    
    ~Solver()
    {
    }
};

template <class T>
class ForwardEuler : public Solver<T>
{
public:
    using Solver<T>::Solver;   // For inheriting the constructor of Solver.
    // void advance(arma::vec (*f)(arma::vec u, double t), int k);
    void advance(T object, int k)
    {   /*
        One step with Forward Euler.

        Parameters
        ----------
        object : object of class T
            Used for accessing the right hand side of the ODE as well
            as the position and velocity data.

        k : int
            Current step in the integration.
        */

        // arma::vec a = f(pos.col(k), 0);
        // this->pos.col(k).print();
        // std::cout << std::endl;
        arma::vec a = object.acceleration(this->pos.col(k), 0);

        this->pos.col(k+1) = this->pos.col(k) + this->dt*this->vel.col(k);
        this->vel.col(k+1) = this->vel.col(k) + this->dt*a;    
    }

};

template <class T>
class VelocityVerlet : public Solver<T>
{
public:
    using Solver<T>::Solver;   // For inheriting the constructor of Solver.
    void advance(T object, int k)
    {   /*
        One step with Velocity Verlet.

        Parameters
        ----------
        object : object of class T
            Used for accessing the right hand side of the ODE as well as
            the position and velocity data.

        k : int
            Current step in the integration.
        */
        
        arma::vec a1 = object.acceleration(this->pos.col(k), 0);
        // arma::vec a1 = f(pos.col(k), 0);
        this->pos.col(k+1) = this->pos.col(k) + this->dt*this->vel.col(k) + this->dt*this->dt/2*a1;

        arma::vec a2 = object.acceleration(this->pos.col(k+1), 0);
        // arma::vec a2 = f(pos.col(k+1), 0);
        this->vel.col(k+1) = this->vel.col(k) + this->dt/2*(a2 + a1);

    }
};

#endif