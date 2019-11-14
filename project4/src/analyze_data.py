import numpy as np
import matplotlib.pyplot as plt


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
    
    print("\nComparing results for a 2x2 spin matrix.")
    print("-"*43)
    print("       analytical   numerical")
    print(f"<|M|>: {M:9.4f} {abs_M_n:9.4f}")
    print(f"<E>:   {E:9.4f} {E_n:9.4f}")
    print(f"Cv:    {Cv:9.4f} {Cv_n:9.4f}")
    print(f"X:     {X:9.4f} {X_n:9.4f}")
    print()


def task_c():
    """
    Cheat sheet:
    Columns are temperatures. Rows are data points.
    data[:, 0]: all data points for temperature 0.
    data[:, 0]: all data points for temperature 0.
    data[5000:, 0]: all data points except the first 5000 for temperature 0.
    All .npy files are calculated with 1e7 MC iterations.
    """


    filename_energy = "/Users/Jon/Desktop/project4/E_convergence_data_20x20.npy"
    filename_magnet = "/Users/Jon/Desktop/project4/E_convergence_data_20x20.npy"

    energy = np.load(filename_energy)
    magnet = np.load(filename_magnet)

    temperatures = energy[0, :]
    temp = 16
    energy = energy[1:, :]
    magnet = magnet[1:, :]

    selection = slice(3000000, -1, 1)
    MC_values = np.arange(1, 1e7+1, 1)
    E_cum_avg = np.cumsum(energy[selection, temp])/np.arange(1, len(energy[selection, temp]) + 1, 1)
    M_cum_avg = np.cumsum(magnet[selection, temp])/np.arange(1, len(magnet[selection, temp]) + 1, 1)


    plt.plot(MC_values[selection], E_cum_avg, label=f"T: {temperatures[temp]:.1f}")    
    plt.xlabel("MC iterations")
    plt.ylabel("Energy, [?]")
    plt.legend(loc="best")
    plt.show()

    # plt.plot(MC_values[selection], M_cum_avg, label=f"T: {temperatures[temp]:.1f}")    
    # plt.xlabel("MC iterations")
    # plt.ylabel("Magnet, [?]")
    # plt.legend(loc="best")
    # plt.show()


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
    compare_values_task_a_and_b()

    # quick_buizz()
    # quick_hist_buizz()
    # task_c()

    pass