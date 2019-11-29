#include "solver.h"

template <class T>
Solver<T>::Solver(int num_steps_input, int num_stellar_objects_input)
    : pos(3*num_stellar_objects_input, num_steps_input+1),
    vel(3*num_stellar_objects_input, num_steps_input+1)
{
    num_steps = num_steps_input;
    num_stellar_objects = num_stellar_objects_input;
}

template <class T>
void Solver<T>::set_initial_conditions(arma::vec U0)
{

    for (int i = 0; i < 3*num_stellar_objects; i++)
    {
        pos(i, 0) = U0(i);
    }

    for (int i = 0; i < 3*num_stellar_objects; i++)
    {
        vel(i, 0) = U0(i + 3*num_stellar_objects);
    }
}

template <class T>
void Solver<T>::solve(T object, double dt_input)
{
    dt = dt_input;
    double t = 0;  // Currently unused.
    
    for (int k = 0; k < num_steps; k++)
    {
        advance(object, k);
    }
}

template <class T>
void Solver<T>::advance(T object, int k)
{
    // Dummy method.
    std::cout << "NotImplementedError" << std::endl;
}


template <class T>
void Solver<T>::write_to_file()
{
    // defining data files
    std::ofstream tull_mc_tull;
    // tull_mc_tull.open("data_files/tull_mc_tull.txt", std::ios_base::app);
    tull_mc_tull.open("data_files/tull_mc_tull.txt");
    
    for (int i = 0; i < num_steps; i++)
    {
        tull_mc_tull << std::setw(20) << std::setprecision(15);
        tull_mc_tull << dt*i;
        tull_mc_tull << std::setw(30) << std::setprecision(15);
        tull_mc_tull << pos(0, i);
        tull_mc_tull << std::setw(30) << std::setprecision(15);
        tull_mc_tull << pos(1, i);
        tull_mc_tull << std::setw(30) << std::setprecision(15);
        tull_mc_tull << vel(0, i);
        tull_mc_tull << std::setw(30) << std::setprecision(15);
        tull_mc_tull << vel(1, i);
        tull_mc_tull << std::endl;
    }
    tull_mc_tull.close();
}

template <class T>
Solver<T>::~Solver()
{
}


template <class T>
void ForwardEuler<T>::advance(T object, int k)
{   /*
    One step with Forward Euler.

    Parameters
    ----------
    f : function pointer
        The right hand side of the ODE.

    k : int
        Current step in the integration.
    */

    // arma::vec a = f(pos.col(k), 0);
    arma::vec a = object.acceleration(pos.col(k), 0);

    pos.col(k+1) = pos.col(k) + dt*vel.col(k);
    vel.col(k+1) = vel.col(k) + dt*a;    
}

template <class T>
void VelocityVerlet<T>::advance(T object, int k)
{   /*
    One step with Velocity Verlet.

    Parameters
    ----------
    f : function pointer
        The right hand side of the ODE.

    k : int
        Current step in the integration.
    */
    
    arma::vec a1 = object.acceleration(pos.col(k), 0);
    // arma::vec a1 = f(pos.col(k), 0);
    pos.col(k+1) = pos.col(k) + dt*vel.col(k) + dt*dt/2*a1;

    arma::vec a2 = object.acceleration(pos.col(k+1), 0);
    // arma::vec a2 = f(pos.col(k+1), 0);
    vel.col(k+1) = vel.col(k) + dt/2*(a2 + a1);

}