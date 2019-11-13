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


    void mc_iteration_convergence_parallel(double temp,
        int num_temp_iter_energy)
    {   /*
        Run the spin flip a given amount of times. Generate data for
        finding how many iterations is needed for convergence. Keep the
        energy values without averaging.

        Parameters
        ----------
        temp : double
            Temperature of current iteration

        num_temp_iter_energy : int
            Number of temperature that has been calculated times number
            of  Monte Carlo iterations. (Jumps forward to the correct
            number of calculated energies/magnetizations.)
        */

        for (int j = 0; j < mc_iterations; j++)
        {   // loops over n*n spin flips a given amount of times
            // saves relevant data for each iteration

            iterate_spin_flip(temp);
            energy_array[num_temp_iter_energy + j] = total_energy;
            magnet_array[num_temp_iter_energy + j] = total_magnetization;
        }
    }


    void iterate_temperature_convergence_parallel(double initial_temp,
        double final_temp, int temps_per_thread, bool ordered_spins)
    {   /*
        Run the calculation for several temperatures in parallel. Keep
        all energy and magnetization raw values.

        Parameters
        ----------
        initial_temp : double
            Initial temperature.
        
        final_temp : double
            Final temperature.

        temps_per_thread : int
            Number of temperatures every thread will calculate.

        ordered_spins : bool
            For toggling initial ordering of the spins to ordered or
            random.
        */
        
        // Temperature interval and initial temperature for each thread.
        double diff_temp = (final_temp - initial_temp)/(world_size*temps_per_thread);
        double initial_temp_thread = initial_temp + diff_temp*temps_per_thread*world_rank;
        double temp;
 
        energy_array = new double[temps_per_thread*mc_iterations];
        magnet_array = new double[temps_per_thread*mc_iterations];

        int root = 0;   // Main thread.
        double* energy_buffer;  // Buffer for MPI_Gather
        double* magnet_buffer;  // Buffer for MPI_Gather

        if (world_rank == root)
        {   /*
            Only the root thread initializes the buffer arrays to
            save memory.
            */
            energy_buffer = new double[world_size*temps_per_thread*mc_iterations];
            magnet_buffer = new double[world_size*temps_per_thread*mc_iterations];

            std::cout << "mc_iterations: " << mc_iterations
            << ", matrix size: " << n << "x" << n << std::endl << std::endl;
        }

        // Starting main timer.
        std::chrono::steady_clock::time_point t_main_1 = std::chrono::steady_clock::now();

        for (int temp_iteration = 0; temp_iteration < temps_per_thread; temp_iteration++)
        {   // Looping over temperature values.

            if (world_rank == root)
            {
                std::cout << "temperature iteration: " << temp_iteration + 1
                << " of: " << temps_per_thread << std::endl;
            }

            if (ordered_spins)
            {   /*
                Resetting the spin matrix for every temperature to the
                initial ordered state.
                */
                spin.initial_spin(ordered_spins);
            }
            
            else
            {   /*
                Resetting the spin matrix for every temperature to the
                initial random state.
                */
                spin.initial_spin();
            }

            temp = initial_temp_thread + diff_temp*temp_iteration;
            // pre-calculated exponential values
            exp_delta_energy[0]  = std::exp(8*J/temp);
            exp_delta_energy[4]  = std::exp(4*J/temp);
            exp_delta_energy[8]  = 1;
            exp_delta_energy[12] = std::exp(-4*J/temp);
            exp_delta_energy[16] = std::exp(-8*J/temp);

            mc_iteration_convergence_parallel(temp,temp_iteration*mc_iterations);

            if (world_rank == root)
            {   // The root thread prints progress information.
                std::chrono::steady_clock::time_point t_main_2 = std::chrono::steady_clock::now();
                std::chrono::duration<double> main_comp_time  = std::chrono::duration_cast<std::chrono::duration<double> >(t_main_2 - t_main_1);
                
                std::cout << "time since beginning: " << main_comp_time.count()
                << std::endl;
            }
        }

        // Starting gather timer.
        std::chrono::steady_clock::time_point t_gather_1 = std::chrono::steady_clock::now();

        if (world_rank == root)
        {
            std::cout << "\ngathering data from all threads" << std::endl;
        }
 
        // Collects the information from every thread and gathers it into buffer.
        MPI_Gather(energy_array, temps_per_thread*mc_iterations, MPI_DOUBLE,
            energy_buffer, temps_per_thread*mc_iterations, MPI_DOUBLE, root,
            MPI_COMM_WORLD);
        MPI_Gather(magnet_array, temps_per_thread*mc_iterations, MPI_DOUBLE,
            magnet_buffer, temps_per_thread*mc_iterations, MPI_DOUBLE, root,
            MPI_COMM_WORLD);

        if (world_rank == root)
        {
            std::chrono::steady_clock::time_point t_gather_2 = std::chrono::steady_clock::now();
            std::chrono::duration<double> gather_comp_time  = std::chrono::duration_cast<std::chrono::duration<double> >(t_gather_2 - t_gather_1);
            std::cout << "gather completed in: " << gather_comp_time.count() << std::endl;
        }

        // Starting write timer.
        std::chrono::steady_clock::time_point t_write_1 = std::chrono::steady_clock::now();

        if (world_rank == root)
        {   // Root thread writes the data to file.

            std::cout << "\nwriting to file" << std::endl;

            E_convergence_data << "mc_iterations: " << mc_iterations
            << " grid: " << n << std::endl;
            M_convergence_data << "mc_iterations: " << mc_iterations << std::endl;
            
            for (int i = 0; i < temps_per_thread*world_size; i++)
            {   // Header with temperature values.
                E_convergence_data << std::setw(20) << std::setprecision(15);
                E_convergence_data << initial_temp_thread + diff_temp*i;
                M_convergence_data << std::setw(20) << std::setprecision(15);
                M_convergence_data << initial_temp_thread + diff_temp*i;
            }
            
            E_convergence_data << std::endl;
            M_convergence_data << std::endl;

            for (int i = 0; i < mc_iterations; i++)
            {   
                for (int j = 0; j < temps_per_thread*world_size; j++)
                {   /*
                    Some funky indexing to write the data as columns
                    instead of rows. The data are stored in a single
                    long array, [[T1E1, T1E2, ...], [T2E1, T2E2, ...], ...]
                    and the indexing fetches first all the E1 values and
                    writes the data as the first row, then the E2 values
                    as the second row, etc. Each column is then all the
                    energy/magnetization values for each temperature. 
                    */
                   
                    E_convergence_data << std::setw(20) << std::setprecision(15) 
                    << energy_buffer[i+j*mc_iterations];
                    M_convergence_data << std::setw(20) << std::setprecision(15)
                    << magnet_buffer[i+j*mc_iterations];
                }
                
                E_convergence_data << std::endl;
                M_convergence_data << std::endl;
            }
            
            if (world_rank == root)
            {
                std::chrono::steady_clock::time_point t_write_2 = std::chrono::steady_clock::now();
                std::chrono::steady_clock::time_point t_final = std::chrono::steady_clock::now();
                std::chrono::duration<double> write_comp_time  = std::chrono::duration_cast<std::chrono::duration<double> >(t_write_2 - t_write_1);
                std::chrono::duration<double> final_comp_time  = std::chrono::duration_cast<std::chrono::duration<double> >(t_final - t_main_1);
                std::cout << "write completed in: " << write_comp_time.count() << std::endl;
                std::cout << "\ntotal time: " << final_comp_time.count() << std::endl;
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
        int temps_per_thread)
    {   /*
        Run the calculation for several temperatures in parallel.
        Calculate average values on the fly. Do not keep all raw data.

        Parameters
        ----------
        initial_temp : double
            Initial temperature.
        
        final_temp : double
            Final temperature.

        temps_per_thread : int
            Number of temperatures every thread will calculate.
        */
        
        double diff_temp = (final_temp - initial_temp)/(world_size*temps_per_thread);
        double initial_temp_thread = initial_temp + diff_temp*temps_per_thread*world_rank;
        double temp;

        double* sum_total_energy_array = new double[temps_per_thread];
        double* sum_total_energy_squared_array = new double[temps_per_thread];
        double* sum_total_magnetization_array = new double[temps_per_thread];
        double* sum_total_magnetization_absolute_array = new double[temps_per_thread];
        double* sum_total_magnetization_squared_array = new double[temps_per_thread];

        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

        for (int temp_iteration = 0; temp_iteration < temps_per_thread; temp_iteration++)
        {   // looping over temperature values

            temp = initial_temp_thread + diff_temp*temp_iteration;
            // pre-calculated exponential values
            exp_delta_energy[0]  = std::exp(8*J/temp);
            exp_delta_energy[4]  = std::exp(4*J/temp);
            exp_delta_energy[8]  = 1;
            exp_delta_energy[12] = std::exp(-4*J/temp);
            exp_delta_energy[16] = std::exp(-8*J/temp);

            mc_iteration_stable(temp);

            sum_total_energy_array[temp_iteration] = sum_total_energy;
            sum_total_energy_squared_array[temp_iteration] = sum_total_energy_squared;
            sum_total_magnetization_array[temp_iteration] = sum_total_magnetization;
            sum_total_magnetization_absolute_array[temp_iteration] = sum_total_magnetization_absolute;
            sum_total_magnetization_squared_array[temp_iteration] = sum_total_magnetization_squared;
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
                for (int i = 0; i < temps_per_thread; i++)
                {   
                    std::cout << initial_temp_thread << std::endl;
                    ising_model_data << std::setw(20) << std::setprecision(15) << initial_temp_thread + diff_temp*i;
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
    int spin_matrix_dim = 40;
    int mc_iterations = 1e3;
    
    double initial_temp = 1;
    double final_temp = 3;
    double temps_per_thread = 2;

    bool ordered_spins = false;

    time_t seed;
    time(&seed);
    
    // ParallelEnergySolver data_model(spin_matrix_dim, mc_iterations, seed);
    // data_model.iterate_temperature_parallel(initial_temp, final_temp, temps_per_thread);

    ParallelEnergySolver convergence_model(spin_matrix_dim, mc_iterations, seed);
    convergence_model.iterate_temperature_convergence_parallel(initial_temp,
        final_temp, temps_per_thread, ordered_spins);

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