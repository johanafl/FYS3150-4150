import time
import sys  # For testing. Can prob. be removed.
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
mpl.rcParams['agg.path.chunksize'] = 10000  # To plot large amounts of data.


def convert_to_npy():
    """
    Converts a bunch of text files to Numpy binary files
    """

    N = 10
    t_load_1 = time.time()
    print("Loading data...")
    for i in range(N):
        # Loads several runs of both initial state systems.
        E_data_random = np.loadtxt(f"data_files/E_convergence_data_20x20_random_{i}.txt",
            skiprows=2, unpack=False)

        np.save(f"data_files/E_convergence_data_20x20_random_{i}.npy",
            E_data_random)

        E_data_ordered = np.loadtxt(f"data_files/E_convergence_data_20x20_ordered_{i}.txt",
            skiprows=2, unpack=False)

        np.save(f"data_files/E_convergence_data_20x20_ordered_{i}.npy",
            E_data_ordered)

        M_data_random = np.loadtxt(f"data_files/M_convergence_data_20x20_random_{i}.txt",
            skiprows=2, unpack=False)

        np.save(f"data_files/M_convergence_data_20x20_random_{i}.npy",
            M_data_random)

        M_data_ordered = np.loadtxt(f"data_files/M_convergence_data_20x20_ordered_{i}.txt",
            skiprows=2, unpack=False)

        np.save(f"data_files/M_convergence_data_20x20_ordered_{i}.npy",
            M_data_ordered)

        t_load_2 = time.time()
        print(f"Data set {i} loaded and saved.")
        print(f"Time since beginning: {t_load_2 - t_load_1:.3f} seconds.\n")


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


class Task4C:
    """
    Implemented as a class just to group task 4c specifics together.
    """
    def E_and_M_as_a_function_of_MC(self):
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
        t_main_1 = time.time()

        MC_values = np.arange(1, 1e6+1, 1)/1e6  # x values for plot

        E_data_random  = []
        E_data_ordered = []
        M_data_random  = []
        M_data_ordered = []
        E_data_random_avg  = np.zeros((int(1e6), 2))
        E_data_ordered_avg = np.zeros((int(1e6), 2))
        M_data_random_avg  = np.zeros((int(1e6), 2))
        M_data_ordered_avg = np.zeros((int(1e6), 2))
        
        N = 10
        t_load_1 = time.time()
        print("Loading data...")
        for i in range(N):
            # Loads several runs of both initial state systems.
            E_data_random.append(np.load(f"data_files/E_convergence_data_20x20_random_{i}.npy"))
            E_data_random[i] = np.cumsum(E_data_random[i][1:], axis=0)
            E_data_random[i] /= 20*20
            E_data_random[i][:, 0] /= np.arange(1, len(E_data_random[i]) + 1, 1)
            E_data_random[i][:, 1] /= np.arange(1, len(E_data_random[i]) + 1, 1)

            E_data_ordered.append(np.load(f"data_files/E_convergence_data_20x20_ordered_{i}.npy"))
            E_data_ordered[i] = np.cumsum(E_data_ordered[i][1:], axis=0)
            E_data_ordered[i] /= 20*20
            E_data_ordered[i][:, 0] /= np.arange(1, len(E_data_ordered[i]) + 1, 1)
            E_data_ordered[i][:, 1] /= np.arange(1, len(E_data_ordered[i]) + 1, 1)

            M_data_random.append(np.load(f"data_files/M_convergence_data_20x20_random_{i}.npy"))
            M_data_random[i] = np.cumsum(np.abs(M_data_random[i][1:]), axis=0)
            M_data_random[i] /= 20*20
            M_data_random[i][:, 0] /= np.arange(1, len(M_data_random[i]) + 1, 1)
            M_data_random[i][:, 1] /= np.arange(1, len(M_data_random[i]) + 1, 1)

            M_data_ordered.append(np.load(f"data_files/M_convergence_data_20x20_ordered_{i}.npy"))
            M_data_ordered[i] = np.cumsum(np.abs(M_data_ordered[i][1:]), axis=0)
            M_data_ordered[i] /= 20*20
            M_data_ordered[i][:, 0] /= np.arange(1, len(M_data_ordered[i]) + 1, 1)
            M_data_ordered[i][:, 1] /= np.arange(1, len(M_data_ordered[i]) + 1, 1)

            # Summing the random initial state data.
            E_data_random_avg += E_data_random[i]
            M_data_random_avg += M_data_random[i]

            # Summing the ordered initial state data.
            E_data_ordered_avg += E_data_ordered[i]
            M_data_ordered_avg += M_data_ordered[i]

            t_load_2 = time.time()
            print(f"Data set {i+1} of {N} loaded.")
            print(f"Time since beginning: {t_load_2 - t_load_1:.3f} seconds.\n")

        # Averaging the random initial state data.
        E_data_random_avg /= N
        M_data_random_avg /= N

        # Averaging the ordered initial state data.
        E_data_ordered_avg /= N
        M_data_ordered_avg /= N


        # T = 1
        t_fig_1 = time.time()
        print("Generating figure and plotting data...")
        
        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
        fig.text(x=0.455, y=0.035, s=r"MC cycles, $10^6$", fontsize=25)
        fig.text(x=0.47, y=1-0.075, s=r"$\tilde{T} = 1 k_bT/J$", fontsize=25)
        fig.text(x=0.035, y=0.65, s=r"$\tilde{E}, [E/J]$", fontsize=25, rotation="vertical")
        fig.text(x=0.035, y=0.28, s=r"$\tilde{M}, [M/a]$", fontsize=25, rotation="vertical")
        

        # E random.
        # ----------------
        
        # Plotting the averaged data as solid black.
        ax[0, 0].plot(MC_values[::4], E_data_random_avg[::4, 0], color="black")

        for i in range(N):
            # Plotting the individual random data sets with alpha=0.1.
            ax[0, 0].plot(MC_values[::4], E_data_random[i][::4, 0], color="black", linestyle="dashed", alpha=0.1)
        
        ax[0, 0].set_title(r"$\langle \tilde{E} \rangle$ random", fontsize=25)
        ax[0, 0].tick_params(labelsize=25)
        ax[0, 0].set_ylim([-1.99725, -1.99725+0.006])
        ax[0, 0].grid()
        
        # E ordered.
        # ----------------

        # Plotting the averaged data as solid black.
        ax[0, 1].plot(MC_values[::4], E_data_ordered_avg[::4, 0], color="black")

        for i in range(N):
            # Plotting the individual random data sets with alpha=0.1.
            ax[0, 1].plot(MC_values[::4], E_data_ordered[i][::4, 0], color="black", linestyle="dashed", alpha=0.1)
        
        ax[0, 1].set_title(r"$\langle \tilde{E} \rangle$ ordered", fontsize=25)
        ax[0, 1].tick_params(labelsize=25)
        ax[0, 1].set_ylim([-2, -2 + 0.006])
        # ax[0, 1].set_ylim([-1.9982, -1.9982+0.006])
        ax[0, 1].grid()
        
        # M random.
        # ----------------

        # Plotting the averaged data as solid black.
        ax[1, 0].plot(MC_values[::4], M_data_random_avg[::4, 0], color="black")

        for i in range(N):
            # Plotting the individual random data sets with alpha=0.1.
            ax[1, 0].plot(MC_values[::4], M_data_random[i][::4, 0], color="black", linestyle="dashed", alpha=0.1)
        
        ax[1, 0].set_title(r"$\langle |\tilde{M}| \rangle$ random", fontsize=25)
        ax[1, 0].tick_params(labelsize=25)
        ax[1, 0].set_ylim([0.9995 - 0.006, 0.9995])
        ax[1, 0].grid()
        
        # M ordered.
        # ----------------

        # Plotting the averaged data as solid black.
        ax[1, 1].plot(MC_values[::4], M_data_ordered_avg[::4, 0], color="black")

        for i in range(N):
            # Plotting the individual random data sets with alpha=0.1.
            ax[1, 1].plot(MC_values[::4], M_data_ordered[i][::4, 0], color="black", linestyle="dashed", alpha=0.1)
        
        ax[1, 1].set_title(r"$\langle |\tilde{M}| \rangle$ ordered", fontsize=25)
        ax[1, 1].tick_params(labelsize=25)
        ax[1, 1].set_ylim([0.996, 0.996+0.006])
        ax[1, 1].grid()

        t_fig_2 = time.time()
        t_main_2 = time.time()
        print(f"Total time: {t_main_2 - t_main_1:.3f} seconds.\n")
        print(f"Plot T = 1 generated in {t_fig_2 - t_fig_1:.3f} seconds.")
        
        plt.show()

        # T = 2.4

        t_fig_1 = time.time()
        print("Generating figure and plotting data...")
        
        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
        fig.text(x=0.455, y=0.035, s=r"MC cycles, $10^6$", fontsize=25)
        fig.text(x=0.47, y=1-0.075, s=r"$\tilde{T} = 2.4 k_bT/J$", fontsize=25)
        fig.text(x=0.035, y=0.65, s=r"$\tilde{E}, [E/J]$", fontsize=25, rotation="vertical")
        fig.text(x=0.035, y=0.28, s=r"$\tilde{M}, [M/a]$", fontsize=25, rotation="vertical")
        
        # E random.
        # ----------------

        # Plotting the averaged data as solid black.
        ax[0, 0].plot(MC_values[::4], E_data_random_avg[::4, 1], color="black", label=r"$\tilde{T}: 2.4 k_bT$")

        for i in range(N):
            # Plotting the individual random data sets with alpha=0.1.
            ax[0, 0].plot(MC_values[::4], E_data_random[i][::4, 1], color="black", linestyle="dashed", alpha=0.1)
        
        ax[0, 0].set_title(r"$\langle \tilde{E} \rangle$ random", fontsize=25)
        ax[0, 0].tick_params(labelsize=25)
        ax[0, 0].set_ylim([-1.3, -1.3+0.1])
        ax[0, 0].grid()
        
        # E ordered.
        # ----------------

        # Plotting the averaged data as solid black.
        ax[0, 1].plot(MC_values[::4], E_data_ordered_avg[::4, 1], color="black")

        for i in range(N):
            # Plotting the individual ordered data sets with alpha=0.1.
            ax[0, 1].plot(MC_values[::4], E_data_ordered[i][::4, 1], color="black", linestyle="dashed", alpha=0.1)
        
        ax[0, 1].set_title(r"$\langle \tilde{E} \rangle$ ordered", fontsize=25)
        ax[0, 1].tick_params(labelsize=25)
        ax[0, 1].set_ylim([-1.3, -1.3+0.1])
        ax[0, 1].grid()
        
        # M random.
        # ----------------

        # Plotting the averaged data as solid black.
        ax[1, 0].plot(MC_values[::4], M_data_random_avg[::4, 1], color="black")

        for i in range(N):
            # Plotting the individual random data sets with alpha=0.1.
            ax[1, 0].plot(MC_values[::4], M_data_random[i][::4, 1], color="black", linestyle="dashed", alpha=0.1)

        ax[1, 0].set_title(r"$\langle |\tilde{M}| \rangle$ random", fontsize=25)
        ax[1, 0].tick_params(labelsize=25)
        ax[1, 0].set_ylim([0.4, 0.4+0.1])
        ax[1, 0].grid()

        # M ordered.
        # ----------------

        # Plotting the averaged data as solid black.
        ax[1, 1].plot(MC_values[::4], M_data_ordered_avg[::4, 1], color="black")

        for i in range(N):
            # Plotting the individual ordered data sets with alpha=0.1.
            ax[1, 1].plot(MC_values[::4], M_data_ordered[i][::4, 1], color="black", linestyle="dashed", alpha=0.1)

        ax[1, 1].set_title(r"$\langle |\tilde{M}| \rangle$ ordered", fontsize=25)
        ax[1, 1].tick_params(labelsize=25)
        ax[1, 1].set_ylim([0.4, 0.4+0.1])
        ax[1, 1].grid()


        t_fig_2 = time.time()
        print(f"Plot T = 1 generated in {t_fig_2 - t_fig_1:.3f} seconds.")
        
        plt.show()


    def accepted_configurations(self):
        """
        Plot the number of accepted configurations as a function of
        temperature. Display one plot for initial ordered state and one
        for initial random state. Average over 5 runs for each case.
        """

        # Arrays for average values.
        random_config_avg  = np.zeros(24)
        ordered_config_avg = np.zeros(24)
        
        fig0, ax0 = plt.subplots(figsize=(10, 9))
        fig1, ax1 = plt.subplots(figsize=(10, 9))
        axins0 = ax0.inset_axes([0.2, 0.4, 0.4, 0.4])
        axins1 = ax1.inset_axes([0.2, 0.4, 0.4, 0.4])

        N = 5
        for i in range(N):
            # Looping over N files of data.
            random_config, temp = np.loadtxt(f"data_files/E_convergence_data_20x20_random_accepted_configs_{i}.txt",
                skiprows=1, max_rows=2, unpack=False)

            ordered_config, temp = np.loadtxt(f"data_files/E_convergence_data_20x20_ordered_accepted_configs_{i}.txt",
                skiprows=1, max_rows=2, unpack=False)

            # Summing the values for averages.
            random_config_avg += random_config
            ordered_config_avg += ordered_config

            # Plots the non-averaged data and scaling the y axis.
            ax0.plot(temp, random_config/1e7, alpha=0.2, color="black")
            axins0.plot(temp, random_config/1e7, alpha=0.2, color="black")
            ax1.plot(temp, ordered_config/1e7, alpha=0.2, color="black")
            axins1.plot(temp, ordered_config/1e7, alpha=0.2, color="black")
        
        
        # Averaging and scaling the y axis.
        random_config_avg  /= N
        random_config_avg  /= 1e7
        ordered_config_avg /= N
        ordered_config_avg /= 1e7

        # Plots the averaged data.
        ax0.plot(temp, random_config_avg, color="black")
        axins0.plot(temp, random_config_avg, color="black")
        ax1.plot(temp, ordered_config_avg, color="black")
        axins1.plot(temp, ordered_config_avg, color="black")
        
        ax0.set_title(r"Initial random state", fontsize=25)
        ax0.set_xlabel(r"$\tilde{T}, [k_bT/J]$", fontsize=25)
        ax0.set_ylabel(r"Accepted configurations, $10^7$", fontsize=25)
        ax1.set_title(r"Initial ordered state", fontsize=25)
        ax1.set_xlabel(r"$\tilde{T}, [k_bT/J]$", fontsize=25)
        ax1.set_ylabel(r"Accepted configurations, $10^7$", fontsize=25)

        ax0.tick_params(labelsize=25)
        ax0.grid()
        ax1.tick_params(labelsize=25)
        ax1.grid()

        # sub region of the original image
        x1, x2, y1, y2 = 2.20329, 2.2042, 6.09273, 6.0991
        axins0.set_xlim(x1, x2)
        axins0.set_ylim(y1, y2)
        axins0.tick_params(labelsize=19)
        axins0.set_xticklabels(["" ,2.2035, 2.204])
        axins1.set_xlim(x1, x2)
        axins1.set_ylim(y1, y2)
        axins1.tick_params(labelsize=19)
        axins1.set_xticklabels(["" ,2.2035, 2.204])

        ax0.indicate_inset_zoom(axins0)
        ax1.indicate_inset_zoom(axins1)
        plt.show()


    def accepted_configurations_longer_interval(self):
        """
        Quick check of a larger data set with more temperatures.
        Update after checking the data: no longer exponential change
        for temperatures over 2.4, but all runs behave almost identically.
        """
        N = 6
        for i in range(N):
            random_config, temp = np.loadtxt(f"data_files/E_convergence_data_20x20_random_accepted_configs_longer_{i}.txt",
                skiprows=1, max_rows=2, unpack=False)

            plt.plot(temp, random_config)
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


def task_4d():

    filename = "data_files/E_convergence_data_20x20_random_0.npy"
    state = filename[36:43]
    data = np.load(filename)

    temperatures = data[0, :]
    data = data[1:, :]

    t0_convergence = int(5e3)   # Number of MC cycles before convergence.
    t0_selection = slice(t0_convergence, None, 1)
    t1_convergence = int(4e5)   # Number of MC cycles before convergence.
    t1_selection = slice(t1_convergence, None, 1)

    std1 = np.std(data[t0_selection, 0])
    std2 = np.std(data[t1_selection, 1])
    std1_scaled = np.std(data[t0_selection, 0]/(20*20))
    std2_scaled = np.std(data[t1_selection, 1]/(20*20))

    print(f"Initial state: {state}.")
    print("-"*22)
    print("std1: ", std1)
    print("std2: ", std2)

    print("std1: ", std1_scaled)
    print("std2: ", std2_scaled)

    bins = np.arange(-802, 802+1, 4)/(20*20)
    data /= 20*20

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))

    # ----------------------
    # T = 1
    n, _, _ = ax[0, 0].hist(data[t0_selection, 0], color="grey", ec="k", bins=bins)
    ax[0, 0].grid()
    ax[0, 0].tick_params(labelsize=20)
    ax[0, 0].set_title(r"$\tilde{T} = 1 k_bT/J$", fontsize=25)
    ax[0, 0].set_ylabel("occurrence", fontsize=25)

    energy_array = np.arange(-800, 800+1, 4)/(20*20)
    MC = 1e6

    ax[1, 0].plot(energy_array, n/(MC - t0_convergence), color="black")
    ax[1, 0].grid()
    ax[1, 0].tick_params(labelsize=20)
    ax[1, 0].set_xlabel(r"$\tilde{E}/(spins), [E/J]$", fontsize=25)
    ax[1, 0].set_ylabel(r"$P(\tilde{E})$", fontsize=25)

    # ----------------------
    # T = 2
    n, _, _ = ax[0, 1].hist(data[t1_selection, 1], color="grey", ec="k", bins=bins)
    ax[0, 1].grid()
    ax[0, 1].tick_params(labelsize=20)
    ax[0, 1].set_title(r"$\tilde{T} = 2.4 k_bT/J$", fontsize=25)

    energy_array = np.arange(-800, 800+1, 4)/(20*20)
    MC = 1e6

    ax[1, 1].plot(energy_array, n/(MC - t1_convergence), color="black")
    ax[1, 1].grid()
    ax[1, 1].tick_params(labelsize=20)
    ax[1, 1].set_xlabel(r"$\tilde{E}/(spins), [E/J]$", fontsize=25)
    # ax[1, 1].set_ylabel("P(E)", fontsize=25)

    plt.show()


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
    
    q = Task4C()
    q.E_and_M_as_a_function_of_MC()
    # q.accepted_configurations()
    # q.accepted_configurations_longer_interval()

    # task_4d()
    # task_4e()

    # quick_buizz()
    # quicker_buizz()
    # quick_n_dirty()

    pass