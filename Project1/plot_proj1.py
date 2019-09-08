import matplotlib.pyplot as plt
import numpy as np
import sys
import os

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


    



def plot_func(axis, n, args=False, exact=False, filename=False, method=True):

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
        axis[0].set_ylabel("u(x)/v(x)", fontsize=18)
        axis[1].loglog(n, np.max(error), "o", label="exact, n={:d}".format(n))
        axis[1].set_title(r"Relative error", fontsize=18)
        axis[1].set_xlabel("n", fontsize=18)
        axis[1].set_ylabel(r"$\epsilon_{max}(n)$", fontsize=18)

        if exact:
            axis[0].plot(x, u, label="u(x)")
    
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
    num_of_args = len(sys.argv)

    arg = np.array(["Methheds", "Epstasy"])
    run = arg[1]

    if run == arg[0]:
        fig, ax   = plt.subplots(1, 2)
        filenames = np.array(["Poisson_values_n_{:d}.txt",
                              "Poisson_values_spes_n_{:d}.txt",
                              "Poisson_values_LU_n_{:d}.txt"])
        methods   = np.array(["General Thomas n={:d}", "Special Thomas n={:d}",
                              "LU armadillo n={:d}"])

        num_n = 3
        exact = False

        for j in range(3):
            exact=False
            for i in range(1, num_n + 1):
                if i == num_n:
                    exact = True
                plot_func(ax, 10**i, exact=exact, filename=filenames[j],
                            method=methods[j])
    
            ax[0].legend()
            plt.show()
    
    else:
        fig, ax  = plt.subplots()
        
        num_n      = 71
        n_array    = np.linspace(1e1, 1e7, 10)

        error_gen  = np.zeros(num_n-10)
        error_spes = np.zeros(num_n-10)

        bar = IncrementalBar("Processing", max=61)

        errors = open("relative_errors.txt", "w")
        errors.write("Thomas              |   Thomas special\n")

        for i in range(10, num_n):
            str = "./proj1_test.out {:d}".format(int(10**(i/10.0)))
            os.system(str)
            values_gen    = np.loadtxt("Poisson_values_n_{:d}.txt"
                                    .format(int(10**(i/10.0))), skiprows=1)
            values_spes   = np.loadtxt("Poisson_values_spes_n_{:d}.txt"
                                    .format(int(10**(i/10.0))), skiprows=1)
            error_gen[i-10]  = np.max(values_gen[:, 2])
            error_spes[i-10] = np.max(values_spes[:, 2])
            errors.write(str(error_gen[i-10]) + "   " + str(error_spes[i-10])
                         + "\n")

            bar.next()
        bar.finish()
        errors.close()

        relative_errors = np.loadtxt("relative_errors.txt", skiprows=1)
        error_general   = relative_errors[:, 0]
        error_special   = relative_errors[:, 1]
        
        for i in range(2):
            plot_func(ax, n_array, args=error_general, filename=error_special,
                        method=False)
    
            ax.legend(loc="best")
            plt.show()
        #plt.savefig("$Maximum relative error per n")
