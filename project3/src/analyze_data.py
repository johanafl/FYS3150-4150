import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def analyze_data():
    """
    Loads data generated of the program gauss_legendre_quadrature.cpp and plots
    for visualization.
    """

    # reads from file
    N_0, error_0, calculated_0, exact_0, comp_time_0 = \
        np.loadtxt("legendre_data.txt", skiprows=1, unpack=True)

    N_1, error_1, calculated_1, exact_1, comp_time_1 = \
        np.loadtxt("laguerre_data.txt", skiprows=1, unpack=True)


    # N vs error plot
    _, ax = plt.subplots(figsize=(10, 8))

    ax.plot(N_0, error_0, label="legendre")
    ax.plot(N_1, error_1, label="laguerre")
    
    ax.set_title("error", fontsize=20)
    ax.set_xlabel("grid points", fontsize=20)
    ax.set_ylabel("error", fontsize=20)
    
    ax.legend(fontsize=20)
    ax.grid()
    ax.tick_params(labelsize=30)
    
    plt.show()

    # N vs computation time plot
    _, ax = plt.subplots(figsize=(10, 8))

    ax.plot(N_0, comp_time_0, label="legendre")
    ax.plot(N_1, comp_time_1, label="laguerre")
    
    ax.set_title("computation time", fontsize=20)
    ax.set_xlabel("grid points", fontsize=20)
    ax.set_ylabel("computation time [s]", fontsize=20)
    
    ax.legend(fontsize=20)
    ax.grid()
    ax.tick_params(labelsize=30)
    
    plt.show()


def analyze_contour_data():

    ranges = np.loadtxt("legendre_contour_data.txt", skiprows=1, max_rows=1)
    grid = np.loadtxt("legendre_contour_data.txt", skiprows=2)
    grid = np.log10(grid)

    N_range = np.arange(*ranges[0:3])
    lambda_range = np.arange(*ranges[3:])


    X, Y = np.meshgrid(N_range, lambda_range)

    plt.contourf(X, Y, grid)
    plt.xlabel("N")
    plt.ylabel("lambda")
    cbar = plt.colorbar()
    # cbar.set_label(r"$log_{10}$ error", fontsize=40)
    # cbar.ax.tick_params(labelsize=30) 
    plt.show()


def integrand(x1, x2):

    tol = 1e-10
    y1 = y2 = z1 = z2 = 0
    
    r1 = np.sqrt(x1**2 + y1**2 + z1**2)
    r2 = np.sqrt(x2**2 + y2**2 + z2**2)


    
    return np.exp(-2*2*(r1 + r2))/(np.abs(r1 - r2))

def plot_function():

    x1 = np.linspace(-1, 1, 1000)
    x2 = 0#np.linspace(-5, 5, 1000)

    plt.plot(x1, integrand(x1, x2))
    plt.show()

    # X, Y = np.meshgrid(x1, x2)

    # Z = integrand(X, Y)

    # Z[np.where(Z > 10)] = 0



    # fig = plt.figure()
    # ax = fig.gca(projection='3d')

    # ax.plot_surface(X, Y, Z)
    # # ax.set_zlim([0, 0.01])
    # plt.show()



if __name__ == "__main__":
    # analyze_data()
    # analyze_contour_data()
    plot_function()
    pass