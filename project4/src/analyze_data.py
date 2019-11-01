import numpy as np
import matplotlib.pyplot as plt

def analyse_energy_const_temp(energy):
    plt.plot(energy, label="Total energy, [?]")
    plt.xlabel("MC iterations")
    plt.ylabel("Energy, [?]")
    plt.legend(loc="best")
    plt.show()

def analyse_magnet_const_temp(magnet):
    plt.plot(magnet, label="Total magnetization, [?]")
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

if __name__ == "__main__":
    filename_energy = "data_files/E_data.txt"
    filename_magnet = "data_files/M_data.txt"

    energy = np.loadtxt(filename_energy)
    magnet = np.loadtxt(filename_magnet)
    # filename_magnet = "M_data.txt"

    analyse_energy_const_temp(energy)
    analyse_magnet_const_temp(magnet)
    pass