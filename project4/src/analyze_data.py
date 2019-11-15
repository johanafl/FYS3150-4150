import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 10000


def compare_values_task_a_and_b():
    """
    Compares analytical values from task a) with the numerically
    calculated values from b).
    """

    J  = 1
    kb = 1
    T  = 1

    # Unpacking numerical values. All values are mean values except T.
    T, E_n, E_squared_n, M_n, M_squared_n, abs_M_n = \
        np.loadtxt("data_files/ising_model_data_2x2.txt", skiprows=2)

    Cv_n = (E_squared_n - E_n**2)/(kb*T)     # Numerical heat capacity.
    X_n  = (M_squared_n - abs_M_n**2)/(kb*T) # Numerical susceptibility.


    def analytical_mean_absolute_magnetic_moment_2x2(T=1):
        """
        The analytical result for the mean absolute value of the
        magnetic moment for a 2x2 spin matrix.

        Parameters
        ----------
        T : float, np.ndarray
            Temperature in units [kb*T*J**-1].

        Returns
        -------
        M : float, numpy.ndarray
            Mean absolute value of the magnetic moment.
        """

        beta = kb*T

        M  = (4*np.exp(8*J*beta) + 8)
        M /= (np.exp(8*J*beta) + np.exp(-8*J*beta) + 6)

        return M


    def analytical_mean_energy_2x2(T=1):
        """
        The analytical result for the mean energy for a 2x2 spin matrix.

        Parameters
        ----------
        T : float, np.ndarray
            Temperature in units [kb*T*J**-1].

        Returns
        -------
        E : float, numpy.ndarray
            Mean energy.
        """

        beta = kb*T

        E  = (8*np.exp(-8*J*beta) - 8*np.exp(8*J*beta))
        E /= (np.exp(8*J*beta) + np.exp(-8*J*beta) + 6)

        return E

    
    def analytical_specific_heat_capacity_2x2(T=1):
        """
        The analytical result for the specific heat capacity.

        Parameters
        ----------
        T : float, np.ndarray
            Temperature in units [kb*T*J**-1].

        Returns
        -------
        Cv : float, numpy.ndarray
            Specific heat capacity.
        """

        beta = kb*T        

        Cv  = 64*J**2*np.exp(8*J*beta) + 64*J**2*np.exp(-8*J*beta)
        Cv /= np.exp(8*J*beta) + np.exp(-8*J*beta) + 6
        Cv -= analytical_mean_energy_2x2()**2
        Cv *= beta/T

        return Cv

    
    def analytical_susceptibility_2x2(T=1):
        """
        The analytical result for the susceptibility.

        Parameters
        ----------
        T : float, np.ndarray
            Temperature in units [kb*T*J**-1].

        Returns
        -------
        X : float, numpy.ndarray
            Susceptibility.
        """

        beta = kb*T        

        X  = 16*np.exp(8*J*beta) + 16
        X /= np.exp(8*J*beta) + np.exp(-8*J*beta) + 6
        X -= analytical_mean_absolute_magnetic_moment_2x2(T)**2
        X *= beta

        return X

    
    M  = analytical_mean_absolute_magnetic_moment_2x2(T)
    E  = analytical_mean_energy_2x2(T)
    Cv = analytical_specific_heat_capacity_2x2(T)
    X  = analytical_susceptibility_2x2(T)

    M_error  = np.abs(M - abs_M_n)
    E_error  = np.abs(E - E_n)
    Cv_error = np.abs(Cv - Cv_n)
    X_error  = np.abs(X - X_n)
    
    print(f"\nComparing results for a 2x2 spin matrix. MC iterations: 1e9. T: {T}")
    print("-"*67)
    print("            analytical      numerical     error")
    print(f"<|M|>: {M:15.8f} {abs_M_n:15.8f} {M_error:15e}")
    print(f"<E>:   {E:15.8f} {E_n:15.8f} {E_error:15e}")
    print(f"Cv:    {Cv:15.8f} {Cv_n:15.8f} {Cv_error:15e}")
    print(f"X:     {X:15.8f} {X_n:15.8f} {X_error:15e}")
    print()


class TaskC:
    """
    Cheat sheet:
    Columns are temperatures. Rows are data points.
    data[:, 0]: all data points for temperature 0.
    data[5000:, 16]: all data points except the first 5000 for temperature 16.
    All .npy files are calculated with 1e7 MC iterations.

    All plots are contained, with their specific parameters, inside
    functions in this function.
    """

    filename_energy_random = "/Users/Jon/Desktop/project4/E_convergence_data_20x20.npy"
    filename_magnet_random = "/Users/Jon/Desktop/project4/E_convergence_data_20x20.npy"
    filename_energy_ordered = "/Users/Jon/Desktop/project4/E_convergence_data_20x20_ordered.npy"
    filename_magnet_ordered = "/Users/Jon/Desktop/project4/E_convergence_data_20x20_ordered.npy"

    MC_values = np.arange(1, 1e7+1, 1)  # x values for plot
    y_scale = 20*20    # Number of spins.

    E_random_data_loaded = False
    M_random_data_loaded = False
    E_ordered_data_loaded = False
    M_ordered_data_loaded = False

    
    def M_cumulative_average_T_24_ordered(self):
        """
        Specific parameters for showing the cumulative average
        magnetization for T = 2.4 from ordered initial state.

        The y axis is scaled with 1/(n*n) (divided by the number of spins)
        such that the comparison of the burn-in time of the different
        scenarios are more easily compared.
        """

        if not self.M_ordered_data_loaded:
            load_time_1 = time.time()
            # Loads the data if it is not already loaded.
            self.magnet_ordered = np.load(self.filename_magnet_ordered)
            self.temperatures_ordered = self.magnet_ordered[0, :]
            self.magnet_ordered = self.magnet_ordered[1:, :]

            load_time_2 = time.time()
            
            print(f"M ordered data loaded in: {load_time_2 - load_time_1:.3f} seconds.")
            self.M_ordered_data_loaded = True
        
        
        temp = 28   # Index for chosen temperature.
        selection = slice(int(0), None, 1)
        M_cum_avg = np.cumsum(self.magnet_ordered[selection, temp])/\
            np.arange(1, len(self.magnet_ordered[selection, temp]) + 1, 1)

        plot_time_1 = time.time()

        fig, ax = plt.subplots(figsize=(10, 8))
        ax.plot(self.MC_values[selection], M_cum_avg/\
            self.y_scale, label=f"T: {self.temperatures_ordered[temp]:.1f}")
        ax.set_ylim([2.60, 2.65])
        ax.set_xlabel("MC iterations", fontsize=30)
        ax.set_ylabel("Magnetization, [?]", fontsize=30)
        ax.legend(loc="best", fontsize=20)
        ax.tick_params(labelsize=25)
        ax.grid()
        
        plt.tight_layout()

        plot_time_2 = time.time()
        print(f"E random data plotted and shown in: {plot_time_2 - plot_time_1:.3f} seconds.")

        plt.show()
    

    def M_cumulative_average_T_24_random(self):
        """
        Specific parameters for showing the cumulative average
        magnetization for T = 2.4 from random initial state.

        The y axis is scaled with 1/(n*n) (divided by the number of spins)
        such that the comparison of the burn-in time of the different
        scenarios are more easily compared.
        """

        if not self.M_random_data_loaded:
            load_time_1 = time.time()
            # Loads the data if it is not already loaded.
            self.magnet_random = np.load(self.filename_magnet_random)
            self.temperatures_random = self.magnet_random[0, :]
            self.magnet_random = self.magnet_random[1:, :]

            load_time_2 = time.time()
            
            print(f"M random data loaded in: {load_time_2 - load_time_1:.3f} seconds.")
            self.M_random_data_loaded = True

        temp = 17
        selection = slice(int(0), None, 1)
        M_cum_avg = np.cumsum(self.magnet_random[selection, temp])/\
            np.arange(1, len(self.magnet_random[selection, temp]) + 1, 1)

        plot_time_1 = time.time()

        fig, ax = plt.subplots(figsize=(10, 8))
        ax.plot(self.MC_values[selection], M_cum_avg/self.y_scale,
            label=f"T: {self.temperatures_random[temp]:.1f}") 
        ax.set_ylim([-2.0-0.5, -2.0])
        ax.set_xlabel("MC iterations", fontsize=30)
        ax.set_ylabel("Magnetization, [?]", fontsize=30)
        ax.legend(loc="best", fontsize=20)
        ax.tick_params(labelsize=25)
        ax.grid()
        plt.tight_layout()

        plot_time_2 = time.time()
        print(f"E random data plotted and shown in: {plot_time_2 - plot_time_1:.3f} seconds.")
        plt.show()


    def E_cumulative_average_T_24_random(self):
        """
        Specific parameters for showing the cumulative average energy
        for T = 2.4 from random initial state.

        The y axis is scaled with 1/(n*n) (divided by the number of spins)
        such that the comparison of the burn-in time of the different
        scenarios are more easily compared.
        """

        if not self.E_random_data_loaded:
            load_time_1 = time.time()
            # Loads the data if it is not already loaded.
            self.energy_random = np.load(self.filename_energy_random)
            self.temperatures_random = self.energy_random[0, :]
            self.energy_random = self.energy_random[1:, :]

            load_time_2 = time.time()
            
            print(f"E random data loaded in: {load_time_2 - load_time_1:.3f} seconds.")
            self.E_random_data_loaded = True

        temp = 17
        selection = slice(int(0), None, 1)
        E_cum_avg = np.cumsum(self.energy_random[selection, temp])/\
            np.arange(1, len(self.energy_random[selection, temp]) + 1, 1)

        plot_time_1 = time.time()
        fig, ax = plt.subplots(figsize=(10, 8))
        
        ax.plot(self.MC_values[selection], E_cum_avg/\
            self.y_scale, label=f"T: {self.temperatures_random[temp]:.1f}")
        
        ax.set_ylim([-2.5, -2.5+0.5])
        ax.set_xlabel("MC iterations", fontsize=30)
        ax.set_ylabel("Energy, [?]", fontsize=30)
        ax.legend(loc="best", fontsize=20)
        ax.tick_params(labelsize=25)
        ax.grid()
        plt.tight_layout()
        
        plot_time_2 = time.time()
        print(f"E random data plotted and shown in: {plot_time_2 - plot_time_1:.3f} seconds.")
        
        plt.show()


    def E_raw_T_24_random(self):
        """
        Specific parameters for showing the raw energy data T = 2.4
        from random initial state.

        The y axis is scaled with 1/(n*n) (divided by the number of spins)
        such that the comparison of the burn-in time of the different
        scenarios are more easily compared.
        """

        if not self.E_random_data_loaded:
            load_time_1 = time.time()
            # Loads the data if it is not already loaded.
            self.energy_random = np.load(self.filename_energy_random)
            self.temperatures_random = self.energy_random[0, :]
            self.energy_random = self.energy_random[1:, :]

            load_time_2 = time.time()
            
            print(f"E random data loaded in: {load_time_2 - load_time_1:.3f} seconds.")
            self.E_random_data_loaded = True

        temp = 17
        selection = slice(int(0), None, 20000)

        plot_time_1 = time.time()
        fig, ax = plt.subplots(figsize=(10, 8))
        
        ax.plot(self.MC_values[selection], self.energy_random[selection, temp]/\
            self.y_scale, label=f"T: {self.temperatures_random[temp]:.1f}")
        
        ax.set_xlabel("MC iterations", fontsize=30)
        ax.set_ylabel("Energy, [?]", fontsize=30)
        ax.legend(loc="best", fontsize=20)
        ax.tick_params(labelsize=25)
        ax.grid()
        
        plt.tight_layout()

        plot_time_2 = time.time()
        print(f"E random data plotted and shown in: {plot_time_2 - plot_time_1:.3f} seconds.")

        plt.show()




    



def analyse_magnet_const_temp(magnet):
    try: 
        temps = magnet[:,0]
        magnet = magnet[:,1:]
        # magnet /= 20**2
        print(np.min(magnet))

        for i in range(np.shape(magnet)[0]): 
            plt.plot(np.cumsum(magnet[i])/np.arange(1, len(magnet[i]) + 1, 1), label=f"T: {temps[i]}")

    except IndexError:
        temps = magnet[0]
        magnet = magnet[1:]
        # magnet /= 20**2
        magnet = np.abs(magnet)
        print(np.min(magnet))
        plt.plot(np.cumsum(magnet)/np.arange(1, len(magnet) + 1, 1), label=f"T: {temps}")
        # plt.plot(magnet, label=f"T: {temps}")
    
    plt.xlabel("MC iterations")
    plt.ylabel("Magnetization, [?]")
    plt.legend(loc="best")
    plt.show()


def analyse_energy(temp, avg_energy):
    plt.plot(temp, avg_energy, label="avgage energy, [?]")
    plt.xlabel(r"Temperature, [$k_{b}T/J$]")
    plt.ylabel("Energy, [?]")
    plt.legend(loc="best")
    plt.show()


def analyse_magnet(temp, avg_magnet):
    plt.plot(temp, avg_magnet, label="avgage magnetization, [?]")
    plt.xlabel(r"Temperature, [$k_{b}T/J$]")
    plt.ylabel("Magnetization, [?]")
    plt.legend(loc="best")
    plt.show()


def analyse_heat_capacity(temp, avg_energy, avg_energy_square):
    plt.plot(temp, avg_energy_square - avg_energy**2, label="?")
    plt.xlabel(r"Temperature, [$k_{b}T/J$]")
    plt.ylabel("Heat capacity, [?]")
    plt.legend(loc="best")
    plt.show()





def quick_buizz():
    data = np.loadtxt("data_files/ising_model_data.txt", skiprows=2, unpack=True)

    T, E, E_squared, M, M_squared, M_abs = data
    kb = 1
    
    Cv = (E_squared - E**2)/(kb*T**2)     # Numerical heat capacity.
    X  = (M_squared - M_abs**2)/(kb*T) # Numerical susceptibility.

    fig, ax = plt.subplots(ncols=2, nrows=2)

    ax[0, 0].plot(T, E)
    ax[0, 0].set_title("E")

    ax[0, 1].plot(T, M_abs)
    ax[0, 1].set_title("M abs")

    ax[1, 0].plot(T, Cv)
    ax[1, 0].set_title("Cv")

    ax[1, 1].plot(T, X)
    ax[1, 1].set_title("X")



    plt.show()


def quick_hist_buizz():

    filename = "data_files/E_convergence_data.txt"
    with open(filename, "r") as infile:
        MC = float(infile.readline().split()[1])

    data = np.loadtxt(filename, skiprows=1, unpack=False)
    # std1, std2 = np.std(data[:, 5001:], axis=1)
    std1 = np.std(data[5000:, 0])
    std2 = np.std(data[5000:, 1])
    std1_scaled = np.std(data[5000:, 0]/(20*20))
    std2_scaled = np.std(data[5000:, 0]/(20*20))

    print("std1: ", std1)
    print("std2: ", std2)

    print("std1: ", std1_scaled)
    print("std2: ", std2_scaled)

    bins = np.arange(-802, 802+1, 4)/(20*20)

    data /= 20*20

    n, _, _ = plt.hist(data[5000:, 0], bins=bins)
    plt.xlabel("E")
    plt.ylabel("occurrence")
    plt.show()

    energy_array = np.arange(-800, 800+1, 4)

    plt.plot(energy_array, n/(MC - 5000))
    plt.xlabel("E")
    plt.ylabel("P(E)")
    plt.show()



if __name__ == "__main__":
    # compare_values_task_a_and_b()
    q = TaskC()
    # q.M_cumulative_average_T_24_random()
    # q.E_cumulative_average_T_24_random()
    q.E_raw_T_24_random()

    # quick_buizz()
    # quick_hist_buizz()

    pass