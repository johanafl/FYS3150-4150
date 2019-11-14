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

            // mc_iteration_convergence_parallel(temp,temp_iteration*mc_iterations);
            for (int j = 0; j < mc_iterations; j++)
            {   /*
                Loops over n*n spin flips a given amount of times and
                saves relevant data for each iteration.
                */
                iterate_spin_flip(temp);
                energy_array[temp_iteration*mc_iterations + j] = total_energy;
                magnet_array[temp_iteration*mc_iterations + j] = total_magnetization;
            }

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
        int temps_per_thread, bool ordered_spins)
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
        int root = 0;   // Root thread.

        double* sum_total_energy_array = new double[temps_per_thread];
        double* sum_total_energy_squared_array = new double[temps_per_thread];
        double* sum_total_magnetization_array = new double[temps_per_thread];
        double* sum_total_magnetization_absolute_array = new double[temps_per_thread];
        double* sum_total_magnetization_squared_array = new double[temps_per_thread];

        // Starting main timer.
        std::chrono::steady_clock::time_point t_main_1 = std::chrono::steady_clock::now();


        if (world_rank == root)
        {   // Root thread prints progress info.

            std::cout << "mc_iterations: " << mc_iterations
            << ", matrix size: " << n << "x" << n << std::endl << std::endl;
        }

        for (int temp_iteration = 0; temp_iteration < temps_per_thread; temp_iteration++)
        {   // looping over temperature values

            if (world_rank == root)
            {   // Root thread prints progress info.
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

            mc_iteration_stable(temp);

            sum_total_energy_array[temp_iteration] = sum_total_energy;
            sum_total_energy_squared_array[temp_iteration] = sum_total_energy_squared;
            sum_total_magnetization_array[temp_iteration] = sum_total_magnetization;
            sum_total_magnetization_absolute_array[temp_iteration] = sum_total_magnetization_absolute;
            sum_total_magnetization_squared_array[temp_iteration] = sum_total_magnetization_squared;

            if (world_rank == root)
            {   // The root thread prints progress information.
                std::chrono::steady_clock::time_point t_main_2 = std::chrono::steady_clock::now();
                std::chrono::duration<double> main_comp_time  = std::chrono::duration_cast<std::chrono::duration<double> >(t_main_2 - t_main_1);
                
                std::cout << "time since beginning: " << main_comp_time.count()
                << std::endl;
            }
        }

        if (world_rank == root)
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
        
        // Starting write timer.
        std::chrono::steady_clock::time_point t_write_1 = std::chrono::steady_clock::now();
        
        for (int rank = 0; rank < world_size; rank++)
        {   // Writing data to file.
            if (rank == world_rank)
            {   
                for (int i = 0; i < temps_per_thread; i++)
                {   
                    // std::cout << initial_temp_thread << std::endl;
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



            if (world_rank == root)
            {   // Root thread prints progress info.

                std::chrono::steady_clock::time_point t_write_2 = std::chrono::steady_clock::now();
                std::chrono::steady_clock::time_point t_final = std::chrono::steady_clock::now();
                std::chrono::duration<double> write_comp_time  = std::chrono::duration_cast<std::chrono::duration<double> >(t_write_2 - t_write_1);
                std::chrono::duration<double> final_comp_time  = std::chrono::duration_cast<std::chrono::duration<double> >(t_final - t_main_1);
                std::cout << "write completed in: " << write_comp_time.count() << std::endl;
                std::cout << "\ntotal time: " << final_comp_time.count() << std::endl;
            }
        
    }


    void iterate_temperature_parallel_more_data(double initial_temp, double final_temp,
        int nr_temps, bool ordered_spins)
    {   /*
        Run the calculation for several temperatures in parallel.
        Calculate average values on the fly. Do not keep all raw data.

        Parameters
        ----------
        initial_temp : double
            Initial temperature.
        
        final_temp : double
            Final temperature.

        nr_temps : int
            Number of temperatures.
        */
        
        double diff_temp = (final_temp - initial_temp)/nr_temps;
        double temp;

        int root = 0;
        double total_gather_time = 0;

        double* sum_total_energy_array = new double[nr_temps];
        double* sum_total_energy_squared_array = new double[nr_temps];
        double* sum_total_magnetization_array = new double[nr_temps];
        double* sum_total_magnetization_absolute_array = new double[nr_temps];
        double* sum_total_magnetization_squared_array = new double[nr_temps];

        double sum_total_energy_all_threads;
        double sum_total_energy_squared_all_threads;
        double sum_total_magnetization_all_threads;
        double sum_total_magnetization_absolute_all_threads;
        double sum_total_magnetization_squared_all_threads;

        // Starting main timer.
        std::chrono::steady_clock::time_point t_main_1 = std::chrono::steady_clock::now();

        if (world_rank == root)
        {   // Root thread prints progress info.

            std::cout << "mc_iterations: " << mc_iterations
            << ", matrix size: " << n << "x" << n << std::endl << std::endl;
            std::cout << "number of threads: " << world_size << std::endl;
        }

        for (int temp_iteration = 0; temp_iteration < nr_temps; temp_iteration++)
        {   // looping over temperature values

            if (world_rank == root)
            {   // Root thread prints progress info.
                std::cout << "temperature iteration: " << temp_iteration + 1
                << " of: " << nr_temps << std::endl;
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

            sum_total_energy_all_threads = 0;
            sum_total_energy_squared_all_threads = 0;
            sum_total_magnetization_all_threads = 0;
            sum_total_magnetization_absolute_all_threads = 0;
            sum_total_magnetization_squared_all_threads = 0;

            temp = initial_temp + diff_temp*temp_iteration;
            // pre-calculated exponential values
            exp_delta_energy[0]  = std::exp(8*J/temp);
            exp_delta_energy[4]  = std::exp(4*J/temp);
            exp_delta_energy[8]  = 1;
            exp_delta_energy[12] = std::exp(-4*J/temp);
            exp_delta_energy[16] = std::exp(-8*J/temp);

            mc_iteration_stable(temp);

            // Starting gather timer.
            std::chrono::steady_clock::time_point t_gather_1 = std::chrono::steady_clock::now();

            MPI_Reduce(&sum_total_energy, &sum_total_energy_all_threads, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&sum_total_energy_squared, &sum_total_energy_squared_all_threads, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&sum_total_magnetization, &sum_total_magnetization_all_threads, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&sum_total_magnetization_absolute, &sum_total_magnetization_absolute_all_threads, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&sum_total_magnetization_squared, &sum_total_magnetization_squared_all_threads, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            // Ending gather timer.
            std::chrono::steady_clock::time_point t_gather_2 = std::chrono::steady_clock::now();
            std::chrono::duration<double> gather_time = std::chrono::duration_cast<std::chrono::duration<double> >(t_gather_2 - t_gather_1);
            total_gather_time += gather_time.count();


            sum_total_energy_array[temp_iteration] = sum_total_energy_all_threads/world_size;
            sum_total_energy_squared_array[temp_iteration] = sum_total_energy_squared_all_threads/world_size;
            sum_total_magnetization_array[temp_iteration] = sum_total_magnetization_all_threads/world_size;
            sum_total_magnetization_absolute_array[temp_iteration] = sum_total_magnetization_absolute_all_threads/world_size;
            sum_total_magnetization_squared_array[temp_iteration] = sum_total_magnetization_squared_all_threads/world_size;

            if (world_rank == root)
            {   // The root thread prints progress information.
                std::chrono::steady_clock::time_point t_main_2 = std::chrono::steady_clock::now();
                std::chrono::duration<double> main_comp_time  = std::chrono::duration_cast<std::chrono::duration<double> >(t_main_2 - t_main_1);
                
                std::cout << "time since beginning: " << main_comp_time.count()
                << std::endl;
            }
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

            for (int i = 0; i < nr_temps; i++)
            {   
                // std::cout << initial_temp << std::endl;
                ising_model_data << std::setw(20) << std::setprecision(15) << initial_temp + diff_temp*i;
                ising_model_data << std::setw(20) << std::setprecision(15) << sum_total_energy_array[i];
                ising_model_data << std::setw(20) << std::setprecision(15) << sum_total_energy_squared_array[i];
                ising_model_data << std::setw(20) << std::setprecision(15) << sum_total_magnetization_array[i];
                ising_model_data << std::setw(20) << std::setprecision(15) << sum_total_magnetization_squared_array[i];
                ising_model_data << std::setw(20) << std::setprecision(15) << sum_total_magnetization_absolute_array[i];
                ising_model_data << std::endl;
            }

            // ending timer
            std::chrono::steady_clock::time_point t_final = std::chrono::steady_clock::now();
            std::chrono::duration<double> comp_time  = std::chrono::duration_cast<std::chrono::duration<double> >(t_final - t_main_1);

            std::cout << "\ntotal gather time: " << total_gather_time << std::endl;
            std::cout << "total time: " << comp_time.count() << std::endl;
        }
            delete[] sum_total_energy_array;
            delete[] sum_total_energy_squared_array;
            delete[] sum_total_magnetization_array;
            delete[] sum_total_magnetization_absolute_array;
            delete[] sum_total_magnetization_squared_array;
        
    }


    ~ParallelEnergySolver()
    {
        MPI_Finalize();
        // delete[] exp_delta_energy;
    }

};


int main()
{   
    // int spin_matrix_dim = 20;
    int spin_matrix_dim = 100;

    // int mc_iterations = 1e4;
    int mc_iterations = 2e6;

    // int stable_iterations = 1e3;
    int stable_iterations = 1e6;
    
    double initial_temp = 2;
    double final_temp = 2.6;
    // double temps_per_thread = 5;
    double nr_temps = 50;

    bool not_ordered_spins = false;
    bool ordered_spins = true;

    time_t seed;
    time(&seed);
    
    ParallelEnergySolver data_model(spin_matrix_dim, mc_iterations, seed);
    data_model.set_stable_iterations(stable_iterations);
    // data_model.iterate_temperature_parallel(initial_temp, final_temp, temps_per_thread, ordered_spins);
    data_model.iterate_temperature_parallel_more_data(initial_temp, final_temp, nr_temps, not_ordered_spins);

    // ParallelEnergySolver convergence_model(spin_matrix_dim, mc_iterations, seed);
    // convergence_model.iterate_temperature_convergence_parallel(initial_temp,
    //     final_temp, temps_per_thread, ordered_spins);

    return 0;
}
