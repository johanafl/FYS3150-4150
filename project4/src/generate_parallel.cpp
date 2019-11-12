#include "energy_solver.h"
#include <chrono>
#include <mpi.h>


class ParallelEnergySolver: public IsingModel
{
private:
    int world_rank;
    int world_size;

    double* energy_array;
    double* magnet_array;

public:
    ParallelEnergySolver(int spin_mat_dim, int mc_iterations_input, long seed) 
    : IsingModel(spin_mat_dim, mc_iterations_input, seed)
    {
        MPI_Init(NULL, NULL);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    }


    void mc_iteration_convergence_parallel(double temp, int num_temp_iter_energy)
    {   /*
        Runs the spin flip a given amount of times. Generates data for finding
        how many iterations is needed for convergence. Keeps the energy values
        without averaging.

        Parameters
        ----------
        temp : double
            Temperature of current iteration

        num_temp_iter_energy : int
            Number of temperature that has been calculated times number of 
            Monte Carlo iterations. (Jumps forward to the correct number of 
            calculated energies/magnetizations.)
        */

        for (int j = 0; j < mc_iterations; j++)
        {   // loops over n*n spin flips a given amount of times
            // saves relevant data for each iteration

            iterate_spin_flip(temp);
            // energy_array[num_temp_iter_energy + j] = total_energy;
            // magnet_array[num_temp_iter_energy + j] = total_magnetization;

            energy_array[num_temp_iter_energy + j] = j;
            magnet_array[num_temp_iter_energy + j] = j;
        }
    }


    void iterate_temperature_convergence_parallel(double initial_temp, double final_temp, int num_of_temp_divided_by_num_of_threds, bool ordered_spins)
    {   /*
        ????????????????
        */
        
        double diff_temp = (final_temp - initial_temp)/(world_size*num_of_temp_divided_by_num_of_threds);
        double init_temp = initial_temp + diff_temp*num_of_temp_divided_by_num_of_threds*world_rank;
        double fin_temp = initial_temp + diff_temp*num_of_temp_divided_by_num_of_threds*(world_rank + 1);


        energy_array = new double[num_of_temp_divided_by_num_of_threds*mc_iterations];
        magnet_array = new double[num_of_temp_divided_by_num_of_threds*mc_iterations];

        int root = 0;   // Main thread.
        double* energy_buffer;
        double* magnet_buffer;

        if (world_rank == root)
        {
            energy_buffer = new double[world_size*num_of_temp_divided_by_num_of_threds*mc_iterations];
            magnet_buffer = new double[world_size*num_of_temp_divided_by_num_of_threds*mc_iterations];
        }

        double temp;

        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

        for (int num_iteration = 0; num_iteration < num_of_temp_divided_by_num_of_threds; num_iteration++)
        {   // looping over temperature values

            if (ordered_spins)
            {   // Resetting the spin matrix for every temperature.
                spin.initial_spin(ordered_spins);
            }
            
            else
            {
                spin.initial_spin();
            }

            temp = init_temp + diff_temp*num_iteration;
            // pre-calculated exponential values
            exp_delta_energy[0]  = std::exp(8*J/temp);
            exp_delta_energy[4]  = std::exp(4*J/temp);
            exp_delta_energy[8]  = 1;
            exp_delta_energy[12] = std::exp(-4*J/temp);
            exp_delta_energy[16] = std::exp(-8*J/temp);

            mc_iteration_convergence_parallel(temp, num_iteration*mc_iterations);

            if (world_rank == 0)
            {
                std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
                std::chrono::duration<double> comp_time  = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
                std::cout << "rank:" << world_rank << ", time since beginning: " << comp_time.count() << std::endl;
                std::cout << std::endl;
            }
        }

        MPI_Gather(energy_array, num_of_temp_divided_by_num_of_threds*mc_iterations, MPI_DOUBLE, energy_buffer, num_of_temp_divided_by_num_of_threds*mc_iterations, MPI_DOUBLE, root, MPI_COMM_WORLD);
        MPI_Gather(magnet_array, num_of_temp_divided_by_num_of_threds*mc_iterations, MPI_DOUBLE, magnet_buffer, num_of_temp_divided_by_num_of_threds*mc_iterations, MPI_DOUBLE, root, MPI_COMM_WORLD);

        if (world_rank == root)
        {
            for (int i = 0; i < num_of_temp_divided_by_num_of_threds*world_size; i++)
            {
                E_convergence_data << std::setw(20) << std::setprecision(15);
                E_convergence_data << init_temp + diff_temp*i;
                M_convergence_data << std::setw(20) << std::setprecision(15);
                M_convergence_data << init_temp + diff_temp*i;
            }
            
            E_convergence_data << std::endl;
            M_convergence_data << std::endl;

            for (int i = 0; i < mc_iterations; i++)
            {   
                for (int j = 0; j < num_of_temp_divided_by_num_of_threds*world_size; j++)
                {
                    E_convergence_data << std::setw(20) << std::setprecision(15) << energy_buffer[i+j*mc_iterations];
                    M_convergence_data << std::setw(20) << std::setprecision(15) << magnet_buffer[i+j*mc_iterations];
                }
                
                E_convergence_data << std::endl;
                M_convergence_data << std::endl;
            }

        }

        if (world_rank == root)
        {
            delete[] energy_buffer;
            delete[] magnet_buffer;
        }

        delete[] energy_array;
        delete[] magnet_array;
    }


    void iterate_temperature_parallel(double initial_temp, double final_temp,
        int num_of_temp_divided_by_num_of_threds)
    {   /*
        ????????????????
        */
        
        double diff_temp = (final_temp - initial_temp)/(world_size*num_of_temp_divided_by_num_of_threds);
        double init_temp = initial_temp + diff_temp*num_of_temp_divided_by_num_of_threds*world_rank;
        double fin_temp = initial_temp + diff_temp*num_of_temp_divided_by_num_of_threds*(world_rank + 1);

        double* sum_total_energy_array = new double[num_of_temp_divided_by_num_of_threds];
        double* sum_total_energy_squared_array = new double[num_of_temp_divided_by_num_of_threds];
        double* sum_total_magnetization_array = new double[num_of_temp_divided_by_num_of_threds];
        double* sum_total_magnetization_absolute_array = new double[num_of_temp_divided_by_num_of_threds];
        double* sum_total_magnetization_squared_array = new double[num_of_temp_divided_by_num_of_threds];

        double temp;

        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

        for (int num_iteration = 0; num_iteration < num_of_temp_divided_by_num_of_threds; num_iteration++)
        {   // looping over temperature values

            temp = init_temp + diff_temp*num_iteration;
            // pre-calculated exponential values
            exp_delta_energy[0]  = std::exp(8*J/temp);
            exp_delta_energy[4]  = std::exp(4*J/temp);
            exp_delta_energy[8]  = 1;
            exp_delta_energy[12] = std::exp(-4*J/temp);
            exp_delta_energy[16] = std::exp(-8*J/temp);

            mc_iteration_stable(temp);

            sum_total_energy_array[num_iteration] = sum_total_energy;
            sum_total_energy_squared_array[num_iteration] = sum_total_energy_squared;
            sum_total_magnetization_array[num_iteration] = sum_total_magnetization;
            sum_total_magnetization_absolute_array[num_iteration] = sum_total_magnetization_absolute;
            sum_total_magnetization_squared_array[num_iteration] = sum_total_magnetization_squared;
        }

        if (world_rank == 0)
        {
            ising_model_data << "mc_iterations: " << mc_iterations;
            ising_model_data << " spin_matrix_dim: " << n;
            ising_model_data << std::endl;
            ising_model_data << std::setw(20) << "T";
            ising_model_data << std::setw(20) << "<E>";
            ising_model_data << std::setw(20) << "<E**2>";
            ising_model_data << std::setw(20) << "<M>";
            ising_model_data << std::setw(20) << "<M**2>";
            ising_model_data << std::setw(20) << "<|M|>";
            ising_model_data << std::endl;
        }

        for (int rank = 0; rank < world_size; rank++)
        {
            if (rank == world_rank)
            {   
                for (int i = 0; i < num_of_temp_divided_by_num_of_threds; i++)
                {   
                    std::cout << init_temp << std::endl;
                    ising_model_data << std::setw(20) << std::setprecision(15) << init_temp + diff_temp*i;
                    ising_model_data << std::setw(20) << std::setprecision(15) << sum_total_energy_array[i];
                    ising_model_data << std::setw(20) << std::setprecision(15) << sum_total_energy_squared_array[i];
                    ising_model_data << std::setw(20) << std::setprecision(15) << sum_total_magnetization_array[i];
                    ising_model_data << std::setw(20) << std::setprecision(15) << sum_total_magnetization_squared_array[i];
                    ising_model_data << std::setw(20) << std::setprecision(15) << sum_total_magnetization_absolute_array[i];
                    ising_model_data << std::endl;
                }
            }

            MPI_Barrier(MPI_COMM_WORLD);
        }

        delete[] sum_total_energy_array;
        delete[] sum_total_energy_squared_array;
        delete[] sum_total_magnetization_array;
        delete[] sum_total_magnetization_absolute_array;
        delete[] sum_total_magnetization_squared_array;

        // ending timer
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        std::chrono::duration<double> comp_time  = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);

        if (world_rank == 0)
        {
            // std::cout << "iterations: " << mc_iterations;
            std::cout << " time: " << comp_time.count() << std::endl;
        }
        
    }

    ~ParallelEnergySolver()
    {
        MPI_Finalize();
        // delete[] exp_delta_energy;
    }

};


int main()
{   
    int spin_matrix_dim = 2;
    int mc_iterations = 1e4;
    
    double initial_temp = 2;
    double final_temp = 2.4;
    double num_of_temperatures = 2;

    bool ordered_spins = false;

    time_t seed;
    time(&seed);
    
    // ParallelEnergySolver data_model(spin_matrix_dim, mc_iterations, seed);
    // data_model.iterate_temperature_parallel(initial_temp, final_temp, num_of_temperatures);

    ParallelEnergySolver convergence_model(spin_matrix_dim, mc_iterations, seed);
    convergence_model.iterate_temperature_convergence_parallel(initial_temp, final_temp, num_of_temperatures, ordered_spins);

    return 0;
}



/*
NB!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
I have noticed that the spin matrix is not initialized for every temperature, 
and thus perhaps introducing covariance and faster convergence than it should. 
We also need to make it possible to initialize with all spins up for every 
temperature.

Make initial_spin() take bool variable such that spins can be ordered.

On line 95-115: Jon -> you are now free to choose how the data should be 
written. Write temp to the right.
*/