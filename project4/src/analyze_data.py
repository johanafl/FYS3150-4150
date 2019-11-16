import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
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

    [Cv] = kb
    """

    filename_energy_random = "data_files/E_convergence_data_20x20_random.txt"
    filename_magnet_random = "data_files/M_convergence_data_20x20_random.txt"
    filename_energy_ordered = "data_files/E_convergence_data_20x20_ordered.txt"
    filename_magnet_ordered = "data_files/M_convergence_data_20x20_ordered.txt"

    MC_values = np.arange(1, 1e6+1, 1)/1e6  # x values for plot
    y_scale = 20*20    # Number of spins.

    tol = 0.02

    E_random_data_loaded = False
    M_random_data_loaded = False
    E_ordered_data_loaded = False
    M_ordered_data_loaded = False

    
    def M_cumulative_average_ordered(self, T=2.4):
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
            self.magnet_ordered = np.loadtxt(self.filename_magnet_ordered, skiprows=2)
            self.temperatures_ordered = self.magnet_ordered[0, :]
            self.magnet_ordered = self.magnet_ordered[1:, :]

            load_time_2 = time.time()
            
            print(f"M ordered data loaded in: {load_time_2 - load_time_1:.3f} seconds.")
            self.M_ordered_data_loaded = True

        # Index for chosen temperature.
        temp_idx = np.where(np.abs(self.temperatures_ordered - T) < self.tol)[0][0]

        print(self.temperatures_ordered[temp_idx])
        
        selection = slice(int(0), None, 1)
        M_cum_avg = np.cumsum(self.magnet_ordered[selection, temp_idx])/\
            np.arange(1, len(self.magnet_ordered[selection, temp_idx]) + 1, 1)

        plot_time_1 = time.time()

        fig, ax = plt.subplots(figsize=(10, 8))
        ax.plot(self.MC_values[selection], M_cum_avg/\
            self.y_scale, label=f"T: {self.temperatures_ordered[temp_idx]:.1f}")
        
        if T == 2.4:
            pass
            ax.set_ylim([0.120, 0.120+0.04])

        if T == 1:
            ax.set_ylim([0.9982, 0.9982+0.002])
        
        ax.set_xlabel("MC iterations", fontsize=30)
        ax.set_ylabel("Magnetization", fontsize=30)
        ax.legend(loc="best", fontsize=20)
        ax.tick_params(labelsize=25)
        # ax.set_title("ordered")
        ax.grid()
        # ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        
        plt.tight_layout()

        plot_time_2 = time.time()
        print(f"M ordered data plotted in: {plot_time_2 - plot_time_1:.3f} seconds.")

        plt.show()
    

    def M_cumulative_average_random(self, T=2.4):
        """
        Specific parameters for showing the cumulative average
        magnetization for T = 2.4, [T] = kbT/J, from random initial
        state.

        The y axis is scaled with 1/(n*n) (divided by the number of spins)
        such that the comparison of the burn-in time of the different
        scenarios are more easily compared.
        """

        if not self.M_random_data_loaded:
            load_time_1 = time.time()
            # Loads the data if it is not already loaded.
            self.magnet_random = np.loadtxt(self.filename_magnet_random, skiprows=2)
            self.temperatures_random = self.magnet_random[0, :]
            self.magnet_random = self.magnet_random[1:, :]

            load_time_2 = time.time()
            
            print(f"M random data loaded in: {load_time_2 - load_time_1:.3f} seconds.")
            self.M_random_data_loaded = True

        # Index for chosen temperature.
        temp_idx = np.where(np.abs(self.temperatures_random - T) < self.tol)[0][0]

        print(self.temperatures_random[temp_idx])
        
        selection = slice(int(0), None, 1)
        M_cum_avg = np.cumsum(self.magnet_random[selection, temp_idx])/\
            np.arange(1, len(self.magnet_random[selection, temp_idx]) + 1, 1)

        plot_time_1 = time.time()

        fig, ax = plt.subplots(figsize=(10, 8))
        ax.plot(self.MC_values[selection], M_cum_avg/\
            self.y_scale, label=f"T: {self.temperatures_random[temp_idx]:.1f}") 
        
        if T == 2.4:
            pass
            # ax.set_ylim([0.08, 0.28])

        if T == 1:
            pass
            ax.set_ylim([0.998, 1])
        
        ax.set_xlabel("MC iterations", fontsize=30)
        ax.set_ylabel("Magnetization", fontsize=30)
        ax.legend(loc="best", fontsize=20)
        ax.tick_params(labelsize=25)
        ax.grid()
        # ax.set_title("random")
        plt.tight_layout()

        plot_time_2 = time.time()
        print(f"M random data plotted in: {plot_time_2 - plot_time_1:.3f} seconds.")
        plt.show()


    def E_cumulative_average_ordered(self, T=2.4):
        """
        Specific parameters for showing the cumulative average energy
        for T = 2.4 from ordered initial state.

        The y axis is scaled with 1/(n*n) (divided by the number of spins)
        such that the comparison of the burn-in time of the different
        scenarios are more easily compared.
        """

        if not self.E_ordered_data_loaded:
            load_time_1 = time.time()
            # Loads the data if it is not already loaded.
            self.energy_ordered = np.loadtxt(self.filename_energy_ordered, skiprows=2)
            self.temperatures_ordered = self.energy_ordered[0, :]
            self.energy_ordered = self.energy_ordered[1:, :]

            load_time_2 = time.time()
            
            print(f"E ordered data loaded in: {load_time_2 - load_time_1:.3f} seconds.")
            self.E_ordered_data_loaded = True

        # Index for chosen temperature.
        temp_idx = np.where(np.abs(self.temperatures_ordered - T) < self.tol)[0][0]

        selection = slice(int(0), None, 1)
        E_cum_avg = np.cumsum(self.energy_ordered[selection, temp_idx])/\
            np.arange(1, len(self.energy_ordered[selection, temp_idx]) + 1, 1)

        plot_time_1 = time.time()
        fig, ax = plt.subplots(figsize=(10, 8))
        
        ax.plot(self.MC_values[selection], E_cum_avg/\
            self.y_scale, label=f"T: {self.temperatures_ordered[temp_idx]:.1f}")
        
        if T == 2.4:
            pass
            # ax.set_ylim([2.5, 2.7])

        if T == 1:
            pass
            ax.set_ylim([-1.998, -1.998+0.002])
        
        ax.set_xlabel("MC iterations", fontsize=30)
        # ax.set_ylabel("Energy, [E/J]", fontsize=30)
        ax.set_ylabel("Energy", fontsize=30)
        ax.legend(loc="best", fontsize=20)
        ax.tick_params(labelsize=25)
        ax.grid()
        # ax.set_title("ordered")
        plt.tight_layout()
        
        plot_time_2 = time.time()
        print(f"E ordered data plotted in: {plot_time_2 - plot_time_1:.3f} seconds.")
        
        plt.show()


    def E_cumulative_average_random(self, T=2.4):
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
            self.energy_random = np.loadtxt(self.filename_energy_random, skiprows=2)
            self.temperatures_random = self.energy_random[0, :]
            self.energy_random = self.energy_random[1:, :]

            load_time_2 = time.time()
            
            print(f"E random data loaded in: {load_time_2 - load_time_1:.3f} seconds.")
            self.E_random_data_loaded = True

        # Index for chosen temperature.
        temp_idx = np.where(np.abs(self.temperatures_random - T) < self.tol)[0][0]

        selection = slice(int(0), None, 1)
        E_cum_avg = np.cumsum(self.energy_random[selection, temp_idx])/\
            np.arange(1, len(self.energy_random[selection, temp_idx]) + 1, 1)

        plot_time_1 = time.time()
        fig, ax = plt.subplots(figsize=(10, 8))
        
        ax.plot(self.MC_values[selection], E_cum_avg/\
            self.y_scale, label=f"T: {self.temperatures_random[temp_idx]:.1f}")
        
        if T == 2.4:
            pass
            # ax.set_ylim([-2.4, -2.2])

        if T == 1:
            pass
            ax.set_ylim([-1.99725, -1.99725+0.002])
        
        ax.set_xlabel("MC iterations", fontsize=30)
        ax.set_ylabel("Energy", fontsize=30)
        ax.legend(loc="best", fontsize=20)
        ax.tick_params(labelsize=25)
        ax.grid()
        # ax.set_title("random")
        plt.tight_layout()
        
        plot_time_2 = time.time()
        print(f"E random data plotted in: {plot_time_2 - plot_time_1:.3f} seconds.")
        
        plt.show()


    def E_raw_random(self, T=2.4):
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
            self.energy_random = np.loadtxt(self.filename_energy_random, skiprows=2)
            self.temperatures_random = self.energy_random[0, :]
            self.energy_random = self.energy_random[1:, :]

            load_time_2 = time.time()
            
            print(f"E random data loaded in: {load_time_2 - load_time_1:.3f} seconds.")
            self.E_random_data_loaded = True

        # Index for chosen temperature.
        temp_idx = np.where(np.abs(self.temperatures_random - T) < self.tol)[0][0]
        
        selection = slice(int(0), None, 20000)

        plot_time_1 = time.time()
        fig, ax = plt.subplots(figsize=(10, 8))
        
        ax.plot(self.MC_values[selection], self.energy_random[selection, temp_idx]/\
            self.y_scale, label=f"T: {self.temperatures_random[temp_idx]:.1f}")
        
        ax.set_xlabel("MC iterations", fontsize=30)
        ax.set_ylabel("Energy, [E/J]", fontsize=30)
        ax.legend(loc="best", fontsize=20)
        ax.tick_params(labelsize=25)
        ax.set_title("random")
        ax.grid()
        
        plt.tight_layout()

        plot_time_2 = time.time()
        print(f"E random data plotted in: {plot_time_2 - plot_time_1:.3f} seconds.")

        plt.show()


    def check_averages(self):
        if not self.E_random_data_loaded:
            load_time_1 = time.time()
            # Loads the data if it is not already loaded.
            self.energy_random = np.loadtxt(self.filename_energy_random, skiprows=2)
            self.temperatures_random = self.energy_random[0, :]
            self.energy_random = self.energy_random[1:, :]

            load_time_2 = time.time()
            
            print(f"E random data loaded in: {load_time_2 - load_time_1:.3f} seconds.")
            self.E_random_data_loaded = True

        

        plt.plot(self.temperatures_random, np.mean(self.energy_random[600000:, :], axis=0))
        plt.legend()
        plt.title("<E>")
        plt.show()

        Cv = (np.mean(self.energy_random[600000:, :]**2, axis=0) - np.mean(self.energy_random[600000:, :], axis=0)**2)/self.temperatures_random**2
        plt.plot(self.temperatures_random, Cv, label="Cv")
        plt.legend()
        plt.title("Cv")
        plt.show()

        if not self.M_random_data_loaded:
            load_time_1 = time.time()
            # Loads the data if it is not already loaded.
            self.magnet_random = np.load(self.filename_magnet_random)
            self.temperatures_random = self.magnet_random[0, :]
            self.magnet_random = self.magnet_random[1:, :]

            load_time_2 = time.time()
            
            print(f"M random data loaded in: {load_time_2 - load_time_1:.3f} seconds.")
            self.M_random_data_loaded = True

        
        plt.plot(self.temperatures_random, np.mean(np.abs(self.magnet_random[600000:, :]), axis=0))
        plt.legend()
        plt.title("<M>")
        plt.show()
        
        X = (np.mean(self.magnet_random[600000:, :]**2, axis=0) - np.mean(np.abs(self.magnet_random[600000:, :]), axis=0)**2)/self.temperatures_random
        plt.plot(self.temperatures_random, X, label="X")
        plt.legend()
        plt.show()




def analyse_heat_capacity(temp, avg_energy, avg_energy_square):
    plt.plot(temp, avg_energy_square - avg_energy**2, label="?")
    plt.xlabel(r"Temperature, [$k_{b}T/J$]")
    plt.ylabel("Heat capacity, [?]")
    plt.legend(loc="best")
    plt.show()





def quick_buizz():
    data = np.loadtxt("data_files/ising_model_data_11x11.txt", skiprows=2, unpack=True)

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

    filename = "/Users/Jon/Desktop/project4/E_convergence_data_20x20.npy"
    # filename = "data_files/E_convergence_data.npy"
    data = np.load(filename)
    # filename = "data_files/E_convergence_data.txt"
    # with open(filename, "r") as infile:
    #     MC = float(infile.readline().split()[1])

    # data = np.loadtxt(filename, skiprows=1, unpack=False)
    temperatures = data[0, :]
    data = data[1:, :]

    temp = 11
    print(temperatures[temp])

    selection = slice(int(1e6), None, 1)
    # std1, std2 = np.std(data[:, 5001:], axis=1)
    # std1 = np.std(data[selection, 0])
    # std2 = np.std(data[selection, temp])
    # std1_scaled = np.std(data[selection, 0]/(20*20))
    # std2_scaled = np.std(data[selection, temp]/(20*20))


    # print("std1: ", std1)
    # print("std2: ", std2)

    # print("std1: ", std1_scaled)
    # print("std2: ", std2_scaled)

    bins = np.arange(-802, 802+1, 4)/(20*20)
    # bins = np.arange(-4e4-2, 4e4+2+1, 4)/(100*100)

    data /= 20*20
    # data /= 100*100

    n, _, _ = plt.hist(data[selection, temp])#, bins=bins)
    plt.xlabel("E")
    plt.ylabel("occurrence")
    plt.show()

    # energy_array = np.arange(-800, 800+1, 4)

    # plt.plot(energy_array, n/(MC - 5000))
    # plt.xlabel("E")
    # plt.ylabel("P(E)")
    # plt.show()


def quicker_buizz():
    MC_values = np.arange(1, 1e6+1, 1)  # x values for plot

    filename_magnet_ordered = "data_files/M_convergence_data_10x10_ordered.npy"
    filename_magnet_random = "data_files/M_convergence_data_10x10_random.npy"

    magnet_ordered = np.load(filename_magnet_ordered)
    magnet_ordered = magnet_ordered[1:]

    
    M_cum_avg = np.cumsum(magnet_ordered)/np.arange(1, len(magnet_ordered) + 1, 1)


    fig, ax = plt.subplots(figsize=(10, 8))
    ax.plot(MC_values, M_cum_avg)
    
    
    ax.set_xlabel("MC iterations", fontsize=30)
    ax.set_ylabel("Magnetization", fontsize=30)
    ax.legend(loc="best", fontsize=20)
    ax.tick_params(labelsize=25)

    plt.show()

    magnet_random = np.load(filename_magnet_random)
    magnet_random = magnet_random[1:]

    
    M_cum_avg_random = np.cumsum(magnet_random)/np.arange(1, len(magnet_random) + 1, 1)


    fig, ax = plt.subplots(figsize=(10, 8))
    ax.plot(MC_values, M_cum_avg_random)
    
    
    ax.set_xlabel("MC iterations", fontsize=30)
    ax.set_ylabel("Magnetization", fontsize=30)
    ax.legend(loc="best", fontsize=20)
    ax.tick_params(labelsize=25)

    plt.show()

def task_4e():
    data_20 = np.loadtxt("data_files/ising_model_data_20x20.txt", skiprows=2, unpack=True)
    data_40 = np.loadtxt("data_files/ising_model_data_40x40.txt", skiprows=2, unpack=True)
    data_60 = np.loadtxt("data_files/ising_model_data_60x60.txt", skiprows=2, unpack=True)
    data_80 = np.loadtxt("data_files/ising_model_data_80x80.txt", skiprows=2, unpack=True)
    data_100 = np.loadtxt("data_files/ising_model_data_100x100.txt", skiprows=2, unpack=True)


    T_20, E_20, E_squared_20, M_20, M_squared_20, M_abs_20 = data_20
    T_40, E_40, E_squared_40, M_40, M_squared_40, M_abs_40 = data_40
    T_60, E_60, E_squared_60, M_60, M_squared_60, M_abs_60 = data_60
    T_80, E_80, E_squared_80, M_80, M_squared_80, M_abs_80 = data_80
    T_100, E_100, E_squared_100, M_100, M_squared_100, M_abs_100 = data_100
    kb = 1

    Cv_20 = (E_squared_20 - E_20**2)/(kb*T_20**2)     # Numerical heat capacity.
    X_20  = (M_squared_20 - M_abs_20**2)/(kb*T_20) # Numerical susceptibility.

    Cv_40 = (E_squared_40 - E_40**2)/(kb*T_40**2)     # Numerical heat capacity.
    X_40  = (M_squared_40 - M_abs_40**2)/(kb*T_40) # Numerical susceptibility.

    Cv_60 = (E_squared_60 - E_60**2)/(kb*T_60**2)     # Numerical heat capacity.
    X_60  = (M_squared_60 - M_abs_60**2)/(kb*T_60) # Numerical susceptibility.

    Cv_80 = (E_squared_80 - E_80**2)/(kb*T_80**2)     # Numerical heat capacity.
    X_80  = (M_squared_80 - M_abs_80**2)/(kb*T_80) # Numerical susceptibility.

    Cv_100 = (E_squared_100 - E_100**2)/(kb*T_100**2)     # Numerical heat capacity.
    X_100  = (M_squared_100 - M_abs_100**2)/(kb*T_100) # Numerical susceptibility.

    fig, ax = plt.subplots(ncols=2, nrows=2)
    fig.text(x=0.455, y=0.035, s=r"$Temperature, [k_bT]$", fontsize=25)

    ax[0, 0].plot(T_20, E_20/(20*20), "--.", label="20x20")
    ax[0, 0].plot(T_40, E_40/(40*40), "--.", label="40x40")
    ax[0, 0].plot(T_60, E_60/(60*60), "--.", label="60x60")
    ax[0, 0].plot(T_80, E_80/(80*80), "--.", label="80x80")
    ax[0, 0].plot(T_100, E_100/(100*100), "--.", label="100x100")
    ax[0, 0].set_title(r"$\langle E \rangle$", fontsize=25)
    ax[0, 0].legend(fontsize=20)
    ax[0, 0].tick_params(labelsize=25)

    ax[0, 1].plot(T_20, M_abs_20/(20*20), "--.", label="20x20")
    ax[0, 1].plot(T_40, M_abs_40/(40*40), "--.", label="40x40")
    ax[0, 1].plot(T_60, M_abs_60/(60*60), "--.", label="60x60")
    ax[0, 1].plot(T_80, M_abs_80/(80*80), "--.", label="80x80")
    ax[0, 1].plot(T_100, M_abs_100/(100*100), "--.", label="100x100")
    ax[0, 1].set_title(r"$\langle|M|\rangle$", fontsize=25)
    ax[0, 1].tick_params(labelsize=25)

    ax[1, 0].plot(T_20, Cv_20/(20*20), "--.", label="20x20")
    ax[1, 0].plot(T_40, Cv_40/(40*40), "--.", label="40x40")
    ax[1, 0].plot(T_60, Cv_60/(60*60), "--.", label="60x60")
    ax[1, 0].plot(T_80, Cv_80/(80*80), "--.", label="80x80")
    ax[1, 0].plot(T_100, Cv_100/(100*100), "--.", label="100x100")
    ax[1, 0].set_title(r"$C_v$", fontsize=25)
    ax[1, 0].tick_params(labelsize=25)

    ax[1, 1].plot(T_20, X_20/(20*20), "--.", label="20x20")
    ax[1, 1].plot(T_40, X_40/(40*40), "--.", label="40x40")
    ax[1, 1].plot(T_60, X_60/(60*60), "--.", label="60x60")
    ax[1, 1].plot(T_80, X_80/(80*80), "--.", label="80x80")
    ax[1, 1].plot(T_100, X_100/(100*100), "--.", label="100x100")
    ax[1, 1].set_title(r"$\chi$", fontsize=25)
    ax[1, 1].tick_params(labelsize=25)

    plt.show()



if __name__ == "__main__":
    # compare_values_task_a_and_b()
    q = TaskC()
    q.M_cumulative_average_ordered(T=2.4)
    # q.M_cumulative_average_random(T=2.4)
    # q.E_cumulative_average_ordered(T=2.4)
    # q.E_cumulative_average_random(T=2.4)
    # q.E_raw_random(T=1)
    # q.check_averages()

    # quick_buizz()
    # task_4e()
    # quick_hist_buizz()

    # quicker_buizz()

    pass