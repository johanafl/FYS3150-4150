import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import sys
import os

def compare_times():
    """
    Function for plotting and comparing timing values for the different
    algorithms. Reads data from text file computed by accompanying C++ program.
    Computes mean values for timings of same grid size and same algorithm, and
    presents them as a plot.
    """

    timing_data = np.loadtxt("compare_times.txt", skiprows=1)

    runs = int(timing_data[0][0])            # number of runs for each grid size value
    num_grid_values = int(timing_data[1][0]) # number of grid size values
    timing_data = timing_data[2:]            # slicing away 'runs' and 'num_grid_values'

    thomas      = np.zeros(num_grid_values)
    thomas_s    = np.zeros(num_grid_values)
    LU          = np.zeros(num_grid_values)
    mean_vals   = np.zeros(num_grid_values)
    grid_values = np.zeros(num_grid_values)

    for i in range(num_grid_values):
        # loops over the timing data and calculates mean values
        tmp = timing_data[i*(runs+1):(i+1)*(runs+1), :]
        grid_values[i] = int(tmp[0][0])
        thomas[i], thomas_s[i], LU[i] = np.mean(tmp[1:], axis=0)

    LU[np.where(LU == -1)] = np.nan     # all -1 values are values not computed
                                        # set to nan so that plot is unaffected

    slope0, intercept0, _, _, _ = stats.linregress(grid_values, thomas)
    slope1, intercept1, _, _, _ = stats.linregress(grid_values, thomas_s)
    print(intercept0 - intercept1)

    plt.loglog(grid_values, thomas,   "-o", label="Thomas")
    plt.loglog(grid_values, thomas_s, "-o", label="Thomas special")
    plt.loglog(grid_values, LU,       "-o", label="LU")
    
    plt.ylabel("Seconds")
    plt.xlabel("Grid points")
    plt.title("Calculation time, mean of 10 runs")
    
    plt.legend(loc="best")
    plt.grid()
    plt.tight_layout(pad=2)
    
    plt.show()

def visualize_error():
    """
    Function for plotting and comparing error values for the different
    algorithms. Reads data from text file computed by accompanying C++ program.
    """

    filenames = ["thomas_algorithm_error.txt",
                 "thomas_algorithm_special_error.txt", "LU_error.txt"]

    special_values = np.array([10, 100, 1000])
    num_special_values = len(special_values)


    for name in filenames:
        # unpacking error data and number of grid points
             
        fig, ax = plt.subplots(figsize=(10,8))
        n, max_rel_error = np.loadtxt(name, unpack=True)
        
        for i in range(num_special_values):
            # finding the closest grid value to a set of 'special' values we
            # wish to highlight in the plot
            idx = np.argmin(np.abs(n - special_values[i]))
            ax.plot(n[idx], max_rel_error[idx], "ro")
        
        ax.loglog(n, max_rel_error)
        ax.set_xlabel("number of grid points, n", fontsize=18)
        ax.set_ylabel("max relative error", fontsize=18)
        ax.set_title(name[:-10], fontsize=18)
        ax.grid()
        
        plt.tight_layout()
        plt.show()




def visualize_data():
    """
    Function for reading data files, organizing them into arrays, and passing
    the data to the plot function plot_data. The input text files contains
    three columns. Column 1 is exact values, column 2 is computed values, and
    column 3 is relative error.
    """

    filenames   = ["thomas_algorithm_n_", "thomas_algorithm_special_n_", "LU_n_"]
    grid_values = ["10", "100", "1000"]

    for name in filenames:
        # iterating over different algorithms
        fig, ax = plt.subplots(figsize=(10,8))
        
        for value in grid_values:
            # iterating over different number of grid points
            filename = name + value + ".txt"    # fname for specific method and grid points
            data     = np.loadtxt(filename, skiprows=1)
            computed = data[:, 1]

            if value == "1000":
                # using only the best resolution for exact value
                exact = data[:, 0]
                plot_data(ax, computed, exact)
            
            else:
                plot_data(ax, computed)

        ax.legend(fontsize=18)
        plt.tight_layout(pad=2.5)

        ax.set_title(name[:-3], fontsize=18)
        ax.set_xlabel("x", fontsize=18)
        ax.set_ylabel("function value", fontsize=18)
        ax.grid()

        plt.show()

        

def plot_data(ax, computed, exact=False):
    """
    A simple function for automation of plotting of computed and exact values.

    Parameters
    ----------
    ax : matplotlib.axes._subplots.AxesSubplot object
        Object on which to draw plots.
    
    computed : numpy.ndarray
        Array containing computed values for plotting.

    exact : bool, numpy.ndarray
        Array containing exact values to plot. Defaults to False for no plot of
        exact values.
    """
    
    n = len(computed)
    x = np.linspace(0, 1, n)
    ax.plot(x, computed, label=n-2)
    
    if not isinstance(exact, bool):
        ax.plot(x, exact, label="exact")





    



def plot_func(axis, n, args=False, exact=False, filename=False, method=True):
    """
    her skjer noe

    parameters
    ----------
    axis : ax object from plt.subplots()

    n : int
        Number of grid points

    kwargs
    ------
    args : False
        kjører på vanlig måte, hvis ikke er det argumenter fra komandolinjen

    exact : False
        Hvis true plottes ekstakte løsningen for den gitte n verdien

    filename : False
        Filnavnet data skal hentes fra

    method : True
        hvis true plottes metoden 
        hvis false plottes kun relativ error
    """

    try:
        if not args:
            values = np.loadtxt(filename.format(n), skiprows=1)
    except ValueError:
        values = args
    
    if method:

        x     = np.linspace(0, 1, n+2)
        u     = values[:, 0]
        v     = values[:, 1]
        error = values[:, 2]

        axis[0].plot(x, v, "-", label=method.format(n))
        axis[0].set_title("Numerical/exact solution", fontsize=18)
        axis[0].set_xlabel("x", fontsize=18)
        axis[0].set_ylabel("function value", fontsize=18)
        axis[1].loglog(n, np.max(error), "o", label="exact, n={:d}".format(n))
        axis[1].set_title(r"Relative error", fontsize=18)
        axis[1].set_xlabel("number of grid points, n", fontsize=18)
        axis[1].set_ylabel(r"$\epsilon_{max}(n)$", fontsize=18)

        if exact:
            axis[0].plot(x, u, label="exact")
    
    else:
        epsilon_general = args
        epsilon_special = filename

        axis.set_title("Relative error for the general Thomas algorithm")
        axis.semilogy(n, epsilon_general, label=r"Max($\epsilon(n)$)")
        axis.set_xlabel("n")
        axis.set_ylabel("Relative error")
        # axis[1].set_title("Relative error for the special Thomas algorithm")
        # axis[1].semilogy(n, epsilon_special, label=r"Max($\epsilon(n)$)")
        # axis[1].set_xlabel("n")
        # axis[1].set_ylabel("Relative error")

if __name__ == "__main__":
    visualize_data()
    # visualize_error()
    # compare_times()
    pass