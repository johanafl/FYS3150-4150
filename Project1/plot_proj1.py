import matplotlib.pyplot as plt
import numpy as np
import sys

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
    num_of_args = len(sys.argv)
    fig, ax  = plt.subplots(1, 2)
    
    if num_of_args >= 2:
        num_n    = num_of_args - 1
        exact    = False

        for i in range(1, num_n + 1):
            if i == num_n:
                exact = True

            filename = str(sys.argv[i])
            values   = np.loadtxt(filename, skiprows=1)
            n        = len(values[:, 0]) - 2

            plot_func(ax, n, args=values, exact=exact)

    else:
        num_n = 3
        exact = False

        for i in range(1, num_n + 1):
            if i == num_n:
                exact = True
            plot_func(ax, 10**i, exact=exact)
    
    ax[0].legend()
    plt.show()