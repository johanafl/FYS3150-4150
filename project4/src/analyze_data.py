import numpy as np
import matplotlib.pyplot as plt

def analyse_energy_const_temp(energy):
    plt.plot(energy[1], label="Total energy, [?]")
    plt.xlabel("MC iterations")
    plt.ylabel("Energy, [?]")
    plt.legend(loc="best")
    plt.show()

def analyse_magnet_const_temp(magnet):
    plt.plot(magnet[2], label="Total magnetization, [?]")
    plt.xlabel("MC iterations")
    plt.ylabel("Magnetization, [?]")
    plt.legend(loc="best")
    plt.show()

def analyse_energy(temp, avr_energy):
    plt.plot(temp, avr_energy, label="Avrage energy, [?]")
    plt.xlabel(r"Temperature, [$k_{b}T/J$]")
    plt.ylabel("Energy, [?]")
    plt.legend(loc="best")
    plt.show()

def analyse_magnet(temp, avr_magnet):
    plt.plot(temp, avr_magnet, label="Avrage magnetization, [?]")
    plt.xlabel(r"Temperature, [$k_{b}T/J$]")
    plt.ylabel("Magnetization, [?]")
    plt.legend(loc="best")
    plt.show()

def analyse_heat_capasity(temp, avr_energy, avr_energy_square):
    plt.plot(temp, avr_energy_square - avr_energy**2, label="?")
    plt.xlabel(r"Temperature, [$k_{b}T/J$]")
    plt.ylabel("Heat capasity, [?]")
    plt.legend(loc="best")
    plt.show()

def analyse_magnet(temp, avr_magnet, avr_magnet_square):
    plt.plot(temp, avr_magnet_square - avr_magnet**2, label="?")
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
        5th, 6th or 7th place in the list in order to see plot of mean energy,
        mean energy^2, mean magnetization, mean magnetization^2, mean
        absolute magnetization, heat capacity or susceptibility respectively.

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
        heat_capasity  = (E_squared - E**2)/scaling_factor

        ax.plot(T, heat_capasity)
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

if __name__ == "__main__":
    filename_energy = "data_files/E_convergence_data.txt"
    filename_magnet = "data_files/M_convergence_data.txt"
    filename_props  = "data_files/ising_model_data.txt"

    energy = np.loadtxt(filename_energy)
    magnet = np.loadtxt(filename_magnet)
    # filename_magnet = "M_data.txt"

    analyse_energy_const_temp(energy)
    analyse_magnet_const_temp(magnet)

    # fig, ax = plt.subplots()

    # analayz_tampratar(ax, filename_props, [True, False, False, False, False, False, False])
    # plt.show()
    pass