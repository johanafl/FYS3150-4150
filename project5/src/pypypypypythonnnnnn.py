import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    data = np.loadtxt("data_files/solver_data.txt", unpack=True)

    planets = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn",
        "Uranus", "Neptune", "Pluto"]
    for i in range(9):
        plt.plot(data[1 + i*6], data[2 + i*6], label=planets[i])
    
    # plt.plot(data[1], data[2], label="mercury")
    # plt.plot(data[1], data[2], label="venus")
    # plt.plot(data[1], data[2], label="earth")
    # plt.plot(data[1], data[2], label="mars")
    # plt.plot(data[1], data[2], label="jupiter")
    # plt.plot(data[1], data[2], label="saturn")
    # plt.plot(data[1], data[2], label="uranus")
    # plt.plot(data[1], data[2], label="neptune")
    # plt.plot(data[1], data[2], label="pluto")

    # plt.plot(data[1], data[2], label="earth")
    # plt.plot(data[7], data[8], label="jupiter")
    # plt.plot(data[13], data[14], label="uranus")
    plt.legend(loc="best")
    plt.axis("equal")
    plt.show()
