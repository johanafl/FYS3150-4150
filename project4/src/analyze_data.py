import numpy as np
import matplotlib.pyplot as plt

def analayzan(ax, file, mean_E=False, mean_M=False, spec_heat=False, susceptibility=False):
    """
    Visualizes the mean energy, E, mean absolute magnetization, |M|,
    the specific heat, Cv, and the susceptibility, X.

    Parameters
    ----------
    ax : 

    file : 
    """
    E  = file[0]
    M  = np.abs(file[1])
    MC = file[some_index]

    if mean_E:
        ax.set_title("Energy of crystal with zomezhings")
        ax.plot(MC, E, label=r"$\langle E \rangle$")
        ax.set_ylabel("Energy")
    
    if mean_M:
        ax.set_title("Magnetization of crystal with zomezhings")
        ax.plot(MC, M, label=r"$\langle |M| \rangle$")
        ax.set_ylabel("Magnetization")


if __name == "__main__":
    fig, ax = plt.subplots()

    analayzan(ax, COMING SOON!)