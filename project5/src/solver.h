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
    int num_steps, num_stellar_objects;
    double dt;

    arma::mat vel, pos;

    Solver(int num_steps_input, int num_stellar_objects_input) 
    : pos(3*num_stellar_objects_input, num_steps_input+1),
    vel(3*num_stellar_objects_input, num_steps_input+1)
    {
        num_steps = num_steps_input;
        num_stellar_objects = num_stellar_objects_input;
    }
    
    void set_initial_conditions(arma::vec U0)
    {   /*
        Load initial values from input vector into position and
        velocity matrices.

        Parameters
        ----------
        U0 : arma::vec
            A vector of initial conditions.
        */

        for (int i = 0; i < 3*num_stellar_objects; i++)
        {
            pos(i, 0) = U0(i);
        }

        for (int i = 0; i < 3*num_stellar_objects; i++)
        {
            vel(i, 0) = U0(i + 3*num_stellar_objects);
        }
    }
    
    void solve(T object, double dt_input)
    {   /*
        Solve the ODE by looping the appropriate advance method.

        Parameters
        ----------
        object : object of class T
            Object containing the RHS of the ODE.

        dt_input : double
            Time step length.
        */
        
        dt = dt_input;
        
        for (int k = 0; k < num_steps; k++)
        {
            advance(object, k);
        }

    }
    
    virtual void advance(T object, int k)
    {   /*
        Dummy method. This method will be overwritten by the child
        classes.

        Parameters
        ----------
        object : object of class T
            Object containing the RHS of the ODE.

        k : int
            Current step in the integration.
        */
        
        std::cout << "NotImplementedError" << std::endl;
    }

    void write_to_file(std::string filepath)
    {   /*
        Write generated data to file.

        Parameters
        ----------
        filepath : std::string
            Path (and filename) to file.

        Note
        ----
        rows: time, x, y, z, vx, vy, vz, ... (for every planet).
        columns: rows for each time step.
        */
        
        std::ofstream outfile;
        outfile.open(filepath);

        for (int i = 0; i < num_steps; i++)
        {
            outfile << std::setw(20) << std::setprecision(15);
            outfile << dt*i;

            for (int j = 0; j < num_stellar_objects; j++)
            {
                outfile << std::setw(30) << std::setprecision(15);
                outfile << pos(3*j + 0, i);
                outfile << std::setw(30) << std::setprecision(15);
                outfile << pos(3*j + 1, i);
                outfile << std::setw(30) << std::setprecision(15);
                outfile << pos(3*j + 2, i);

                outfile << std::setw(30) << std::setprecision(15);
                outfile << vel(3*j + 0, i);
                outfile << std::setw(30) << std::setprecision(15);
                outfile << vel(3*j + 1, i);
                outfile << std::setw(30) << std::setprecision(15);
                outfile << vel(3*j + 2, i);
            }
            outfile << std::endl;
        }
        outfile.close();
    }
};

template <class T>
class ForwardEuler : public Solver<T>
{
public:
    using Solver<T>::Solver;   // For inheriting the constructor of Solver.
    void advance(T object, int k)
    {   /*
        One step with Forward Euler.

        Parameters
        ----------
        object : object of class T
            Used for accessing the right hand side of the ODE.

        k : int
            Current step in the integration.
        */

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
            Used for accessing the right hand side of the ODE.

        k : int
            Current step in the integration.
        */
        
        arma::vec a1 = object.acceleration(this->pos.col(k), 0);
        this->pos.col(k+1) = this->pos.col(k) + this->dt*this->vel.col(k) + this->dt*this->dt/2*a1;

        arma::vec a2 = object.acceleration(this->pos.col(k+1), 0);
        this->vel.col(k+1) = this->vel.col(k) + this->dt/2*(a2 + a1);

    }
};

#endif