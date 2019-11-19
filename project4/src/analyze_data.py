import time
import sys  # For testing. Can prob. be removed.
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from scipy.interpolate import UnivariateSpline
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


class Task4E:
    def __init__(self):
        """
        Plot <E>, <|M|>, Cv, X for grids 2x2, 20x20, 40x40, 60x60, 80x80,
        100x100.
        """
        data_2   = np.loadtxt("data_files/ising_model_data_2x2.txt", skiprows=2, unpack=True)
        data_20  = np.loadtxt("data_files/ising_model_data_20x20.txt", skiprows=2, unpack=True)
        data_40  = np.loadtxt("data_files/ising_model_data_40x40.txt", skiprows=2, unpack=True)
        data_60  = np.loadtxt("data_files/ising_model_data_60x60.txt", skiprows=2, unpack=True)
        data_80  = np.loadtxt("data_files/ising_model_data_80x80.txt", skiprows=2, unpack=True)
        data_100 = np.loadtxt("data_files/ising_model_data_100x100.txt", skiprows=2, unpack=True)

        self.T_2, self.E_2, self.E_squared_2, self.M_2, self.M_squared_2, self.M_abs_2 = data_2
        self.T_20, self.E_20, self.E_squared_20, self.M_20, self.M_squared_20, self.M_abs_20 = data_20
        self.T_40, self.E_40, self.E_squared_40, self.M_40, self.M_squared_40, self.M_abs_40 = data_40
        self.T_60, self.E_60, self.E_squared_60, self.M_60, self.M_squared_60, self.M_abs_60 = data_60
        self.T_80, self.E_80, self.E_squared_80, self.M_80, self.M_squared_80, self.M_abs_80 = data_80
        self.T_100, self.E_100, self.E_squared_100, self.M_100, self.M_squared_100, self.M_abs_100 = data_100
        
        kb = 1
        self.kb = kb
        self.J = 1

        self.Cv_2 = (self.E_squared_2 - self.E_2**2)/(kb*self.T_2**2)     # Numerical heat capacity.
        self.X_abs_2  = (self.M_squared_2 - self.M_abs_2**2)/(kb*self.T_2) # Numerical susceptibility.
        self.X_2  = (self.M_squared_2 - self.M_2**2)/(kb*self.T_2) # Numerical susceptibility.

        self.Cv_20 = (self.E_squared_20 - self.E_20**2)/(kb*self.T_20**2)     # Numerical heat capacity.
        self.X_abs_20  = (self.M_squared_20 - self.M_abs_20**2)/(kb*self.T_20) # Numerical susceptibility.
        self.X_20  = (self.M_squared_20 - self.M_20**2)/(kb*self.T_20) # Numerical susceptibility.

        self.Cv_40 = (self.E_squared_40 - self.E_40**2)/(kb*self.T_40**2)     # Numerical heat capacity.
        self.X_abs_40  = (self.M_squared_40 - self.M_abs_40**2)/(kb*self.T_40) # Numerical susceptibility.
        self.X_40  = (self.M_squared_40 - self.M_40**2)/(kb*self.T_40) # Numerical susceptibility.

        self.Cv_60 = (self.E_squared_60 - self.E_60**2)/(kb*self.T_60**2)     # Numerical heat capacity.
        self.X_abs_60  = (self.M_squared_60 - self.M_abs_60**2)/(kb*self.T_60) # Numerical susceptibility.
        self.X_60  = (self.M_squared_60 - self.M_60**2)/(kb*self.T_60) # Numerical susceptibility.

        self.Cv_80 = (self.E_squared_80 - self.E_80**2)/(kb*self.T_80**2)     # Numerical heat capacity.
        self.X_abs_80  = (self.M_squared_80 - self.M_abs_80**2)/(kb*self.T_80) # Numerical susceptibility.
        self.X_80  = (self.M_squared_80 - self.M_80**2)/(kb*self.T_80) # Numerical susceptibility.

        self.Cv_100 = (self.E_squared_100 - self.E_100**2)/(kb*self.T_100**2)     # Numerical heat capacity.
        self.X_abs_100  = (self.M_squared_100 - self.M_abs_100**2)/(kb*self.T_100) # Numerical susceptibility.
        self.X_100  = (self.M_squared_100 - self.M_100**2)/(kb*self.T_100) # Numerical susceptibility.

        self.Cv_2_spl   = UnivariateSpline(self.T_2, self.Cv_2)
        self.Cv_20_spl  = UnivariateSpline(self.T_20, self.Cv_20)
        self.Cv_40_spl  = UnivariateSpline(self.T_40, self.Cv_40)
        self.Cv_60_spl  = UnivariateSpline(self.T_60, self.Cv_60)
        self.Cv_80_spl  = UnivariateSpline(self.T_80, self.Cv_80)
        self.Cv_100_spl = UnivariateSpline(self.T_100, self.Cv_100)
        
        # self.Tc_2   = self.T_2[np.argmax(self.Cv_2)]
        # self.Tc_20  = self.T_20[np.argmax(self.Cv_20)]
        # self.Tc_40  = self.T_40[np.argmax(self.Cv_40)]
        # self.Tc_60  = self.T_60[np.argmax(self.Cv_60)]
        # self.Tc_80  = self.T_80[np.argmax(self.Cv_80)]
        # self.Tc_100 = self.T_100[np.argmax(self.Cv_100)]

    def visualize_data(self, show_plot=True, show_spline_plot=False):
        
        fine_temp = np.linspace(2, 2.6, 10000)

        if show_spline_plot:

            fig, ax = plt.subplots(figsize=(10, 8))

            line, = ax.plot(fine_temp, self.Cv_2_spl(fine_temp)/(2*2), label="2x2")
            ax.plot(self.T_2, self.Cv_2/(2*2), ".", label="2x2", color=line.get_color())
            
            line, = ax.plot(fine_temp, self.Cv_20_spl(fine_temp)/(20*20), label="20x20")
            ax.plot(self.T_20, self.Cv_20/(20*20), ".", label="20x20", color=line.get_color())
            
            line, = ax.plot(fine_temp, self.Cv_40_spl(fine_temp)/(40*40), label="40x40")
            ax.plot(self.T_40, self.Cv_40/(40*40), ".", label="40x40", color=line.get_color())
            
            line, = ax.plot(fine_temp, self.Cv_60_spl(fine_temp)/(60*60), label="60x60")
            ax.plot(self.T_60, self.Cv_60/(60*60), ".", label="60x60", color=line.get_color())
            
            line, = ax.plot(fine_temp, self.Cv_80_spl(fine_temp)/(80*80), label="80x80")
            ax.plot(self.T_80, self.Cv_80/(80*80), ".", label="80x80", color=line.get_color())
            
            line, = ax.plot(fine_temp, self.Cv_100_spl(fine_temp)/(100*100), label="100x100")
            ax.plot(self.T_100, self.Cv_100/(100*100), ".", label="100x100", color=line.get_color())

            ax.set_title(r"$\tilde{C}_v$", fontsize=25)
            ax.set_ylabel(r"$\tilde{C}_v/(spins), [k_B\sigma^2_{\tilde{E}}/ \tilde{T}^{2}]$", fontsize=25)
            ax.tick_params(labelsize=25)
            ax.grid()

            plt.show()


        if show_plot:

            fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(10, 8))
            fig.text(x=0.48, y=0.035, s=r"$\tilde{T}, [k_BT/J]$", fontsize=25)

            ax[0, 0].plot(self.T_2, self.E_2/(2*2), "--.", label="2x2")
            ax[0, 0].plot(self.T_20, self.E_20/(20*20), "--.", label="20x20")
            ax[0, 0].plot(self.T_40, self.E_40/(40*40), "--.", label="40x40")
            ax[0, 0].plot(self.T_60, self.E_60/(60*60), "--.", label="60x60")
            ax[0, 0].plot(self.T_80, self.E_80/(80*80), "--.", label="80x80")
            ax[0, 0].plot(self.T_100, self.E_100/(100*100), "--.", label="100x100")
            ax[0, 0].set_title(r"$\langle \tilde{E} \rangle$", fontsize=25)
            ax[0, 0].set_ylabel(r"$\tilde{E}/(spins), [E/J]$", fontsize=25)
            ax[0, 0].legend(fontsize=20)
            ax[0, 0].tick_params(labelsize=25)
            ax[0, 0].grid()

            ax[0, 1].plot(self.T_2, self.M_abs_2/(2*2), "--.", label="2x2")
            ax[0, 1].plot(self.T_20, self.M_abs_20/(20*20), "--.", label="20x20")
            ax[0, 1].plot(self.T_40, self.M_abs_40/(40*40), "--.", label="40x40")
            ax[0, 1].plot(self.T_60, self.M_abs_60/(60*60), "--.", label="60x60")
            ax[0, 1].plot(self.T_80, self.M_abs_80/(80*80), "--.", label="80x80")
            ax[0, 1].plot(self.T_100, self.M_abs_100/(100*100), "--.", label="100x100")
            ax[0, 1].set_title(r"$\langle|\tilde{M}|\rangle$", fontsize=25)
            ax[0, 1].set_ylabel(r"$\tilde{M}/(spins), [M/a]$", fontsize=25)
            ax[0, 1].tick_params(labelsize=25)
            ax[0, 1].grid()

            line, = ax[1, 0].plot(fine_temp, self.Cv_2_spl(fine_temp)/(2*2), label="fitted")
            ax[1, 0].plot(self.T_2, self.Cv_2/(2*2), ".", color=line.get_color(), label="data points")
            
            line, = ax[1, 0].plot(fine_temp, self.Cv_20_spl(fine_temp)/(20*20))
            ax[1, 0].plot(self.T_20, self.Cv_20/(20*20), ".", color=line.get_color())
            
            line, = ax[1, 0].plot(fine_temp, self.Cv_40_spl(fine_temp)/(40*40))
            ax[1, 0].plot(self.T_40, self.Cv_40/(40*40), ".", color=line.get_color())
            
            line, = ax[1, 0].plot(fine_temp, self.Cv_60_spl(fine_temp)/(60*60))
            ax[1, 0].plot(self.T_60, self.Cv_60/(60*60), ".", color=line.get_color())
            
            line, = ax[1, 0].plot(fine_temp, self.Cv_80_spl(fine_temp)/(80*80))
            ax[1, 0].plot(self.T_80, self.Cv_80/(80*80), ".", color=line.get_color())
            
            line, = ax[1, 0].plot(fine_temp, self.Cv_100_spl(fine_temp)/(100*100))
            ax[1, 0].plot(self.T_100, self.Cv_100/(100*100), ".", color=line.get_color())

            ax[1, 0].set_title(r"$\tilde{C}_v$", fontsize=25)
            ax[1, 0].set_ylabel(r"$\tilde{C}_v/(spins), [k_B\sigma^2_{\tilde{E}}/ \tilde{T}^{2}]$", fontsize=25)
            ax[1, 0].tick_params(labelsize=25)
            ax[1, 0].grid()
            ax[1, 0].legend(loc="upper left", fontsize=20)

            ax[1, 1].plot(self.T_2, self.X_abs_2/(2*2), "--.", label="2x2")
            ax[1, 1].plot(self.T_20, self.X_abs_20/(20*20), "--.", label="20x20")
            ax[1, 1].plot(self.T_40, self.X_abs_40/(40*40), "--.", label="40x40")
            ax[1, 1].plot(self.T_60, self.X_abs_60/(60*60), "--.", label="60x60")
            ax[1, 1].plot(self.T_80, self.X_abs_80/(80*80), "--.", label="80x80")
            ax[1, 1].plot(self.T_100, self.X_abs_100/(100*100), "--.", label="100x100")
            ax[1, 1].set_title(r"$\tilde{\chi}_{|\tilde{M}|}$", fontsize=25)
            ax[1, 1].set_ylabel(r"$\tilde{\chi}_{|\tilde{M}|}/(spins), [\sigma^2_{\tilde{|M|}}/(J \tilde{T})]$", fontsize=25)
            ax[1, 1].tick_params(labelsize=25)
            ax[1, 1].grid()

            plt.show()


            # fig, ax = plt.subplots(ncols=1, nrows=2, figsize=(10, 8))
            # fig.text(x=0.48, y=0.035, s=r"$\tilde{T}, [k_bT/J]$", fontsize=25)

            # ax[0].plot(self.T_20, self.M_squared_20/(20*20), "--.", label="20x20")
            # ax[0].plot(self.T_40, self.M_squared_40/(40*40), "--.", label="40x40")
            # ax[0].plot(self.T_60, self.M_squared_60/(60*60), "--.", label="60x60")
            # ax[0].plot(self.T_80, self.M_squared_80/(80*80), "--.", label="80x80")
            # ax[0].plot(self.T_100, self.M_squared_100/(100*100), "--.", label="100x100")
            # ax[0].set_title(r"$\langle \tilde{E} \rangle$", fontsize=25)
            # ax[0].set_ylabel(r"$\tilde{E}/(spins), [E/J]$", fontsize=25)
            # ax[0].legend(fontsize=20)
            # ax[0].tick_params(labelsize=25)
            # ax[0].grid()

            # ax[1].plot(self.T_20, self.X_20/(20*20), "--.", label="20x20")
            # ax[1].plot(self.T_40, self.X_40/(40*40), "--.", label="40x40")
            # ax[1].plot(self.T_60, self.X_60/(60*60), "--.", label="60x60")
            # ax[1].plot(self.T_80, self.X_80/(80*80), "--.", label="80x80")
            # ax[1].plot(self.T_100, self.X_100/(100*100), "--.", label="100x100")
            # ax[1].set_title(r"$\langle|\tilde{M}|\rangle$", fontsize=25)
            # ax[1].set_ylabel(r"$\tilde{M}/(spins), [M/a]$", fontsize=25)
            # ax[1].tick_params(labelsize=25)
            # ax[1].grid()

            # plt.show()

    def approximate_critical_temperature_for_L_infty(self, show_plot=False):
        """
        Analytical solution for the time derivative of Cv for 2x2 to
        calculate the strange constant, a and more.
        """

        # Using the temperature derivative of Cv for 2x2 to find a.
        self.T = np.linspace(1, 3, 10000)
        arg    = 8*self.J/(self.kb*self.T)
        res    = 2*(8*self.J)**3*np.sinh(arg)*(3*np.cosh(arg) + 1)/(self.kb**2*self.T**4*(np.cosh(arg) + 3)**3)
        res   -= 3*(8*self.J)**3*np.sinh(arg)/(self.kb**2*self.T**4*(np.cosh(arg) + 3)**2)
        res   -= 2*(8*self.J)**2*3*(np.cosh(arg) + 1)/(self.kb*self.T**3*(np.cosh(arg) + 3)**2)

        idx = np.argmin(np.abs(res))
        Tc_2x2 = self.T[idx]    # Critical temperature for 2x2.
        Tc = 2*self.J/(np.log(1 + np.sqrt(2))*self.kb)  # Critical temperature for L -> infty.
        self.a = (Tc_2x2 - Tc)/2**(-1)

        fine_temp = np.linspace(2.2, 2.4, 10000)

        # Extracting the indices of the max Cv for each lattice.
        idx_2   = np.argmax(self.Cv_2_spl(fine_temp))
        idx_20  = np.argmax(self.Cv_20_spl(fine_temp))
        idx_40  = np.argmax(self.Cv_40_spl(fine_temp))
        idx_60  = np.argmax(self.Cv_60_spl(fine_temp))
        idx_80  = np.argmax(self.Cv_80_spl(fine_temp))
        idx_100 = np.argmax(self.Cv_100_spl(fine_temp))

        # Finding the critical temperatures for each lattice.
        Tc_2     = fine_temp[idx_2]
        Tc_20    = fine_temp[idx_20]
        Tc_40    = fine_temp[idx_40]
        Tc_60    = fine_temp[idx_60]
        Tc_80    = fine_temp[idx_80]
        Tc_100   = fine_temp[idx_100]
        Tc_exact = 2/np.log(1 + np.sqrt(2))

        print(f"2x2: {Tc_2:.5f}, exact: {Tc_exact:.5f}, diff: {Tc_2 - Tc_exact:.5f}")
        print(f"20x20: {Tc_20:.5f}, exact: {Tc_exact:.5f}, diff: {Tc_20 - Tc_exact:.5f}")
        print(f"40x40: {Tc_40:.5f}, exact: {Tc_exact:.5f}, diff: {Tc_40 - Tc_exact:.5f}")
        print(f"60x60: {Tc_60:.5f}, exact: {Tc_exact:.5f}, diff: {Tc_60 - Tc_exact:.5f}")
        print(f"80x80: {Tc_80:.5f}, exact: {Tc_exact:.5f}, diff: {Tc_80 - Tc_exact:.5f}")
        print(f"100x100: {Tc_100:.5f}, exact: {Tc_exact:.5f}, diff: {Tc_100 - Tc_exact:.5f}")





if __name__ == "__main__":
    # compare_values_task_a_and_b()
    
    # q = Task4C()
    # q.E_and_M_as_a_function_of_MC()
    # q.accepted_configurations()
    # q.accepted_configurations_longer_interval()

    # task_4d()
    q = Task4E()
    # q.visualize_data()
    q.approximate_critical_temperature_for_L_infty()

    # quick_buizz()
    # quicker_buizz()
    # quick_n_dirty()

    pass