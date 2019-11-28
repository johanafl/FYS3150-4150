import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    data = np.loadtxt("data_files/tull_mc_tull.txt", unpack=True)
    plt.plot(data[1],data[2])
    plt.show()