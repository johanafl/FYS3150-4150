import numpy as np
import matplotlib.pyplot as plt

def visualize_eigenvalues():
    """
    Visualizing calculated, exact and error.
    """
    calc, exact, error = np.loadtxt("eigenvalues.txt", skiprows=1, unpack=True)
    N = len(error)

    plt.plot(np.arange(N), error)
    plt.show()


if __name__ == "__main__":
    visualize_eigenvalues()
    pass