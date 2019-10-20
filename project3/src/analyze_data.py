import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

def analyze_leglag_data():
    """
    Loads data generated from the program gauss_legendre_quadrature.cpp and
    gauss_laguerre_quadrature.cpp and plots for visualization.
    """

    path_0 = "data_files/legendre_data.txt"
    path_1 = "data_files/laguerre_data.txt"

    try:
        # reads from file
        N_0, error_0, calculated_0, exact_0, comp_time_0 = \
            np.loadtxt(path_0, skiprows=1, unpack=True)

        N_1, error_1, calculated_1, exact_1, comp_time_1 = \
            np.loadtxt(path_1, skiprows=1, unpack=True)

    except:
        print(f"Files {path_0} and/or {path_1} were not found. Exiting.")
        sys.exit()


    # N vs error plot
    _, ax = plt.subplots(figsize=(10, 8))

    ax.plot(N_0, error_0, label="legendre")
    ax.plot(N_1, error_1, label="laguerre")
    
    ax.set_title("error", fontsize=25)
    ax.set_xlabel("grid points", fontsize=25)
    ax.set_ylabel("error", fontsize=25)
    
    ax.legend(fontsize=20)
    ax.grid()
    ax.tick_params(labelsize=30)
    
    # plt.show()

    # N vs computation time plot
    _, ax = plt.subplots(figsize=(10, 8))

    ax.plot(N_0, comp_time_0, label="legendre")
    ax.plot(N_1, comp_time_1, label="laguerre")
    
    ax.set_title("computation time", fontsize=25)
    ax.set_xlabel("grid points", fontsize=25)
    ax.set_ylabel("computation time [s]", fontsize=25)
    
    ax.legend(fontsize=20)
    ax.grid()
    ax.tick_params(labelsize=30)
    
    # plt.show()
    plt.clf()

    # error vs time plot
    _, ax = plt.subplots(figsize=(10, 8))

    ax.semilogy(comp_time_0, error_0, label="legendre")
    ax.semilogy(comp_time_1, error_1, label="laguerre")
    
    # ax.set_title("lol", fontsize=25)
    ax.set_xlabel("computation time [s]", fontsize=30)
    ax.set_ylabel("error", fontsize=30)
    ax.set_xlim([None, 29])
    # ax.set_ylim([None, 0.05])
    
    ax.legend(fontsize=20)
    ax.grid()
    ax.tick_params(labelsize=30)
    
    plt.show()


def analyze_mc_data():
    """
    Loads data generated of the programs mc_integration.cpp,
    mc_integration_improved.cpp, mc_integration_improved_parallel.cpp and plots
    for visualization.
    """

    path_0 = "data_files/mc_data_avg.txt"
    path_1 = "data_files/mc_data_O3.txt"
    path_2 = "data_files/mc_improved_parallel_O3_data.txt"
    path_3 = "data_files/mc_improved_parallel_data.txt"

    try:
        # reads from file
        N_0, error_0, calculated_0, exact_0, comp_time_0 = \
            np.loadtxt(path_0, skiprows=1, unpack=True)

        N_1, error_1, calculated_1, exact_1, comp_time_1 = \
            np.loadtxt(path_1, skiprows=1, unpack=True)

        N_2, error_2, calculated_2, exact_2, comp_time_2, variance_2 = \
            np.loadtxt(path_2, skiprows=1, unpack=True)

        N_3, error_3, calculated_3, exact_3, comp_time_3, variance_3 = \
            np.loadtxt(path_3, skiprows=1, unpack=True)

    except OSError:
        print(f"Files {path_0} and/or {path_1} and/or {path_2} and/or {path_3}"
        " were not found Exiting.")
        sys.exit()


    # N vs error plot
    _, ax = plt.subplots(figsize=(10, 8))

    ax.plot(N_0, error_0, label="mc avg")
    ax.plot(N_1, error_1, label="mc O3")
    ax.plot(N_2, error_2, label="mc O3 parallel")
    ax.plot(N_3, error_3, label="mc parallel")
    
    ax.set_title("error", fontsize=25)
    ax.set_xlabel("iterations", fontsize=25)
    ax.set_ylabel("error", fontsize=25)
    
    ax.legend(fontsize=20)
    ax.grid()
    ax.tick_params(labelsize=30)
    
    plt.show()

    # N vs computation time plot
    _, ax = plt.subplots(figsize=(10, 8))

    ax.plot(N_0, comp_time_0, label="mc avg")
    ax.plot(N_1, comp_time_1, label="mc O3")
    ax.plot(N_2, comp_time_2, label="mc O3 parallel")
    ax.plot(N_3, comp_time_3, label="mc parallel")
    
    ax.set_title("computation time", fontsize=25)
    ax.set_xlabel("iterations", fontsize=25)
    ax.set_ylabel("computation time [s]", fontsize=25)
    
    ax.legend(fontsize=20)
    ax.grid()
    ax.tick_params(labelsize=30)
    
    plt.show()

    # # error vs time plot
    # _, ax = plt.subplots(figsize=(10, 8))

    # ax.semilogy(comp_time_0, error_0, label="mc avg")
    # ax.semilogy(comp_time_1, error_1, label="mc O3")
    # ax.semilogy(comp_time_2, error_2, label="mc O3 parallel")
    # ax.semilogy(comp_time_3, error_3, label="mc parallel")
    
    # # ax.set_title("lol", fontsize=25)
    # ax.set_xlabel("computation time [s]", fontsize=30)
    # ax.set_ylabel("error", fontsize=30)
    # # ax.set_xlim([None, 29])
    # # ax.set_ylim([None, 0.05])
    
    # ax.legend(fontsize=20)
    # ax.grid()
    # ax.tick_params(labelsize=30)
    
    # plt.show()


def compare_speeds():
    """
    Speed comparison.
    """
    _, ax = plt.subplots(figsize=(10, 8))


    for i in range(2, 10+1, 2):
        path_0 = "data_files/mc_improved_parallel_" + f"{i}" + "thread_data.txt"
        path_1 = "data_files/mc_improved_parallel_" + f"{i}" + "thread_O3_data.txt"

        N_0, error_0, calculated_0, exact_0, comp_time_0, variance_0 = \
            np.loadtxt(path_0, skiprows=1, unpack=True)

        N_1, error_1, calculated_1, exact_1, comp_time_1, variance_1 = \
            np.loadtxt(path_1, skiprows=1, unpack=True)

        ax.plot(N_0, comp_time_0, linestyle="dashed", label=f"{i} threads")
        ax.plot(N_1, comp_time_1, label=f"{i} threads O3")

        ax.set_title("computation time", fontsize=25)
        ax.set_xlabel("iterations", fontsize=25)
        ax.set_ylabel("computation time [s]", fontsize=25)

        # ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
        ax.set_xticklabels(["0", "0", "1e1", "1e2", "1e3", "1e4", "1e5"])
        
        ax.legend(fontsize=10)
        ax.grid()
        ax.tick_params(labelsize=30)

    plt.show()


def analyze_contour_data():
    """
    Reads contour data from file generated by gauss_legendre_quadrature.cpp.
    """

    
    # path = "data_files/legendre_contour_data_high_res.txt"
    path = "data_files/mc_contour_data_high_res.txt"
    
    print(f"Reading file {path}.")

    try:
        ranges = np.loadtxt(path, skiprows=1, max_rows=1)
        grid = np.loadtxt(path, skiprows=2)
    except:
        print(f"The file {path} was not found. Exiting.")
        sys.exit()
    
    grid = np.log10(grid)

    N_start, N_end, dN = ranges[0:3]
    lambda_start, lambda_end, dlambda = ranges[3:]
    
    N_range      = np.arange(N_start, N_end+dN, dN)
    lambda_range = np.arange(lambda_start, lambda_end+dlambda, dlambda)
    X, Y = np.meshgrid(N_range, lambda_range)

    # extracting the indices of the lowest error
    lambda_min, N_min = np.unravel_index(np.argmin(grid), np.shape(grid))
    print(f"Lowest error (log10): {grid[lambda_min, N_min]}")
    print(f"N: {N_range[N_min]}")
    print(f"lambda: {lambda_range[lambda_min]}")


    _, ax = plt.subplots()

    mappable = ax.contourf(X, Y, grid)
    # ax.plot(N_range[N_min], lambda_range[lambda_min], "ro")
    ax.set_title(path)
    ax.set_xlabel("N", fontsize=30)
    ax.set_ylabel("lambda", fontsize=30)
    ax.tick_params(labelsize=30)
    
    cbar = plt.colorbar(mappable)
    cbar.set_label(r"$log_{10}$ error", fontsize=40)
    cbar.ax.tick_params(labelsize=30) 
    
    plt.show()


def integrand(x1, x2):

    y1 = y2 = z1 = z2 = 0
    
    r1 = np.sqrt(x1**2 + y1**2 + z1**2)
    r2 = np.sqrt(x2**2 + y2**2 + z2**2)

    return np.exp(-2*2*(r1 + r2))/np.sqrt( (x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2 )


def plot_integrand():
    """
    Plots the integrand. Two subplots of varying x1 and x2 respectively.
    """

    endpoints = 2
    xlim = 2   # limit for the x axis
    ylim = 1/5
    
    x1 = np.linspace(-endpoints, endpoints, 10000)

    _, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 5))
    
    ax.plot(x1, integrand(x1=x1, x2=0))
    
    ax.set_ylabel("f(x)", fontsize=30)
    ax.set_xlabel("x", fontsize=30)
    ax.set_xlim([-xlim, xlim])
    ax.set_ylim([-0.05, ylim])
    
    ax.tick_params(labelsize=30)
    ax.grid()

    ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    
    plt.tight_layout(pad=2)
    plt.show()

    _, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 5))
    
    ax.plot(x1, integrand(x1=x1, x2=0))
    
    ax.set_ylabel("f(x)", fontsize=30)
    ax.set_xlabel("x", fontsize=30)
    ax.set_xlim([-xlim, xlim])
    
    ax.tick_params(labelsize=30)
    ax.grid()

    ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))

    plt.tight_layout(pad=2)
    plt.show()

def illustrate_distributions(distribution="uniform"):
    """
    Illustrates how the uniform and exponential distributions overlap with the
    integrand.

    Parameters
    ----------
    distribution : str
        String input of either "uniform" or "exponential".
    """

    x0 = np.linspace(-0.01, -1, 1000)
    x1 = np.linspace(0.01, 1, 1000)

    _, ax = plt.subplots()
    
    ax.plot(x0, np.abs(1/x0), color="black", label="f(x)")
    ax.plot(x1, 1/x1, color="black")

    if distribution == "uniform":
        ax.plot([-1, 1], [10, 10], label="uniform distribution")
        ax.plot([-1, -1], [0, 10.7], linestyle="dashed", alpha=0.8, color="black")
        ax.plot([1, 1], [0, 10.7], linestyle="dashed", alpha=0.8, color="black")

    elif distribution == "exponential":
        lam = 15
        ax.plot(x0, 5*lam*np.exp(-lam*x1), color="orange")
        ax.plot(x1, 5*lam*np.exp(-lam*x1), label="exponential distribution", color="orange")
    
    ax.set_ylim([-0, 20])
    ax.set_xlabel("x", fontsize=30)
    ax.set_ylabel("f(x)", fontsize=30)
    ax.grid()
    ax.tick_params(labelsize=30)
    ax.legend(fontsize=15)
    plt.show()

if __name__ == "__main__":
    # analyze_leglag_data()
    # analyze_mc_data()
    # comparke_speeds()
    # analyze_contour_data()
    # plot_integrand()
    illustrate_distributions(distribution="uniform")
    pass