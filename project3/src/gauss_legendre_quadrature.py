import numpy as np
import matplotlib.pyplot as plt

def analyze_data():

    N, error, calculated, exact, comp_time = \
        np.loadtxt("legendre_data.txt", skiprows=1, unpack=True)


    _, ax = plt.subplots(figsize=(10, 8))

    ax.plot(N, error, label="error")
    ax.set_xlabel("grid points", fontsize=20)
    ax.set_ylabel("error", fontsize=20)
    ax.legend(fontsize=20)
    ax.grid()
    ax.tick_params(labelsize=30)
    plt.show()

    _, ax = plt.subplots(figsize=(10, 8))

    ax.plot(N, comp_time, label="comp time")
    ax.set_xlabel("grid points", fontsize=20)
    ax.set_ylabel("computation time [s]", fontsize=20)
    ax.legend(fontsize=20)
    ax.grid()
    ax.tick_params(labelsize=30)
    plt.show()



if __name__ == "__main__":
    analyze_data()
    pass