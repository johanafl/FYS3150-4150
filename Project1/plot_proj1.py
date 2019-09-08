import matplotlib.pyplot as plt
import numpy as np
import sys

def compare_times():
    """
    Function for plotting and comparing timing values for the different
    algorithms.
    """

    data = np.loadtxt("compare_times.txt", skiprows=1)

    runs = int(data[0][0])            # number of runs for each grid size value
    num_grid_values = int(data[1][0]) # number of grid size values
    data = data[2:]                   # slicing away 'runs' and 'num_grid_values'

    thomas      = np.zeros(num_grid_values)
    thomas_s    = np.zeros(num_grid_values)
    LU          = np.zeros(num_grid_values)
    mean_vals   = np.zeros(num_grid_values)
    grid_values = np.zeros(num_grid_values)

    for i in range(num_grid_values):
        tmp = data[i*(runs+1):(i+1)*(runs+1), :]
        grid_values[i] = int(tmp[0][0])
        thomas[i], thomas_s[i], LU[i] = np.mean(tmp[1:], axis=0)

    LU[np.where(LU == -1)] = np.nan     # all -1 values are values not computed

    plt.loglog(grid_values, thomas, label="Thomas")
    plt.loglog(grid_values, thomas_s, label="Thomas special")
    plt.loglog(grid_values, LU, label="LU")
    plt.legend(loc="best")
    plt.ylabel("Seconds")
    plt.xlabel("Grid points")
    plt.grid()
    plt.tight_layout(pad=2)
    plt.show()


    



def plot_func(axis, n, args=False, exact=False):

    try:
        if not args:
            values = np.loadtxt("Poisson_values_n_{:d}.txt".format(n), skiprows=1)
    except ValueError:
        values = args

    x     = np.linspace(0, 1, n+2)
    u     = values[:, 0]
    v     = values[:, 1]
    error = values[:, 2]

    axis[0].plot(x, v, "-", label="v, n={:d}".format(n))
    axis[1].loglog(n, np.max(error), "o", label="exact, n={:d}".format(n))

    if exact:
        axis[0].plot(x, u, label="u(x)")


if __name__ == "__main__":
    compare_times()
    # num_of_args = len(sys.argv)
    # fig, ax  = plt.subplots(1, 2)
    
    # if num_of_args >= 2:
    #     num_n    = num_of_args - 1
    #     exact    = False

    #     for i in range(1, num_n + 1):
    #         if i == num_n:
    #             exact = True

    #         filename = str(sys.argv[i])
    #         values   = np.loadtxt(filename, skiprows=1)
    #         n        = len(values[:, 0]) - 2

    #         plot_func(ax, n, args=values, exact=exact)

    # else:
    #     num_n = 3
    #     exact = False

    #     for i in range(1, num_n + 1):
    #         if i == num_n:
    #             exact = True
    #         plot_func(ax, 10**i, exact=exact)
    
    # ax[0].legend()
    # plt.show()
    pass