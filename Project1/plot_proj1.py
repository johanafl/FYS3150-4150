import matplotlib.pyplot as plt
import numpy as np
import sys

def plot(axis, n, exact=False):
    values = np.loadtxt("Poisson_values_n_{:d}.txt".format(n), skiprows=1)


    x     = np.linspace(0, 1, n+2)
    u     = values[:, 0]
    v     = values[:, 1]
    error = values[:, 2]

    axis[0].plot(x, v, "-", label="v, n={:d}".format(n))
    axis[1].loglog(n, np.max(error), "o", label="{:d}".format(n))

    if exact:
        axis[0].plot(x, u, label="u(x)")


if __name__ == "__main__":
    fig, ax = plt.subplots(1, 2)
    num_n = 3
    exact = False

    for i in range(1, num_n+1):
        if i == num_n:
            exact = True
        plot(ax, 10**i, exact)
    
    ax[0].legend()
    plt.show()