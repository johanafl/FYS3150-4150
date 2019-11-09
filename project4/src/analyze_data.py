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




def analyse_energy_const_temp():

    # filename_energy = "data_files/E_convergence_data_20x20.txt"
    # filename_magnet = "data_files/M_convergence_data_20x20.txt"

    filename_energy = "data_files/E_convergence_data.txt"
    filename_magnet = "data_files/M_convergence_data.txt"

    energy = np.loadtxt(filename_energy, skiprows=1)
    magnet = np.loadtxt(filename_magnet, skiprows=1)

    # plt.plot(np.cumsum(energy[0, :10000])/np.arange(1, len(energy[0, :10000]) + 1, 1))
    plt.plot(np.cumsum(magnet[0])/np.arange(1, len(magnet[0]) + 1, 1))
    # plt.plot(energy[0, :1000])
    
    plt.xlabel("MC iterations")
    plt.ylabel("Magnet, [?]")
    plt.legend(loc="best")
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

def analyse_magnet(temp, avg_magnet, avg_magnet_square):
    plt.plot(temp, avg_magnet_square - avg_magnet**2, label="?")
    plt.xlabel(r"Temperature, [$k_{b}T/J$]")
    plt.ylabel("Suceptibility, [?]")
    plt.legend(loc="best")
    plt.show()

def analayz_tampratar(ax, filename, properties=[]):
    """
    Function for plotting the different properties of the crystal against 
    temperature, i.e., plotting average energy, average (energy^2), ...

    Parameters
    ----------

    ax : axis object

    filename : string
        name of the data file

    properties : bool, list_like(?)
        List of boolian values for triggering the plot for the different
        physical properties; the True value needs to be in 1st, 2nd, 3rd, 4th,
        5th, 6th or 7th place in the properties list in order to see plot of
        mean energy, mean energy^2, mean magnetization, mean magnetization^2,
        mean absolute magnetization, heat capacity or susceptibility
        respectively.

        OBS!!
        Don't send in more than 1 True value!!!
    """
    data = np.loadtxt(filename)

    T         = data[:, 0]
    E         = data[:, 1]
    E_squared = data[:, 2]
    M         = data[:, 3]
    M_squared = data[:, 4]
    M_absolut = data[:, 5]
    
    if properties[0]:
        ax.plot(T, E)
        ax.set_xlabel("Temperature")
        ax.set_ylabel("Energy")
        ax.set_title("Nice title here")
    
    if properties[1]:
        ax.plot(T, E_squared)
        ax.set_xlabel("Temperature")
        ax.set_ylabel("Energy squared")
        ax.set_title("Nice title here")
    
    if properties[2]:
        ax.plot(T, M)
        ax.set_xlabel("Temperature")
        ax.set_ylabel("Magnetization")
        ax.set_title("Nice title here")
    
    if properties[3]:
        ax.plot(T, M_squared)
        ax.set_xlabel("Temperature")
        ax.set_ylabel("Magnetization squared")
        ax.set_title("Nice title here")

    if properties[4]:
        ax.plot(T, M_absolut)
        ax.set_xlabel("Temperature")
        ax.set_ylabel("Absolute magnetization")
        ax.set_title("Nice title here")
    
    if properties[5]:
        scaling_factor = 1
        heat_capacity  = (E_squared - E**2)/scaling_factor

        ax.plot(T, heat_capacity)
        ax.set_xlabel("Temperature")
        ax.set_ylabel("Heat capacity")
        ax.set_title("Nice title here")
    
    if properties[6]:
        scaling_factor = 1
        susceptibility = (M_squared - M**2)/scaling_factor

        ax.plot(T, susceptibility)
        ax.set_xlabel("Temperature")
        ax.set_ylabel("Susceptibility")
        ax.set_title("Nice title here")


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

    data = np.loadtxt(filename, skiprows=1)
    std1, std2 = np.std(data[:, 5001:], axis=1)
    std1_scaled, std2_scaled = np.std(data[:, 5001:]/(20*20), axis=1)

    print("std1: ", std1)
    print("std2: ", std2)

    print("std1: ", std1_scaled)
    print("std2: ", std2_scaled)

    bins = np.arange(-802, 802+1, 4)/(20*20)

    data /= 20*20

    n, _, _ = plt.hist(data[1, 5001:], bins=bins)
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

    # quick_buizz()
    # quick_hist_buizz()

    # # filename_magnet = "M_data.txt"
    analyse_energy_const_temp()
    # analyse_magnet_const_temp(magnet)

    # # fig, ax = plt.subplots()

    # analayz_tampratar(ax, filename_props, [True, False, False, False, False, False, False])
    # plt.show()
    pass