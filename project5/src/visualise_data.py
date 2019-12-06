import numpy as np
import matplotlib.pyplot as plt
import constants as c


def task_5c():
    earth_mass = 5.972e24
    solar_mass   = 1.9891e30
    G = 6.67e-11
    dts = ["0.001000", "0.010000", "0.050000", "0.100000"]
    
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
    
    # Velocity Verlet
    #-----------
    
    filepath = f"data_files/task_5c_vv_dt={dts[0]}.txt"
    data = np.loadtxt(filepath, unpack=True)
    
    tol = float(dts[0])
    idx = np.where( np.abs(data[2]) < tol )
    diff1 = data[1][idx][np.where( (data[1][idx]) > 0 )]

    v = np.sqrt(data[4]**2 + data[5]**2 + data[6]**2)*c.AU*1000/31556926
    r = np.sqrt(data[1]**2 + data[2]**2 + data[3]**2)*c.AU*1000
    K = 1/2*earth_mass*v**2
    V = -G*solar_mass*earth_mass/r
    E1 = K + V

    # print(np.diff(E)[::100]/E[0])

    ax[0, 0].plot(data[1], data[2], label="Earth")
    ax[0, 0].set_title(f"{dts[0]}")
    ax[0, 0].axis("equal")
    
    #-----------

    filepath = f"data_files/task_5c_vv_dt={dts[1]}.txt"
    data = np.loadtxt(filepath, unpack=True)
    
    tol = 3e-2
    idx = np.where( np.abs(data[2]) < tol )
    diff2 = data[1][idx][np.where( (data[1][idx]) > 0 )]

    v = np.sqrt(data[4]**2 + data[5]**2 + data[6]**2)*c.AU*1000/31556926
    r = np.sqrt(data[1]**2 + data[2]**2 + data[3]**2)*c.AU*1000
    K = 1/2*earth_mass*v**2
    V = -G*solar_mass*earth_mass/r
    E2 = K + V

    ax[0, 1].plot(data[1], data[2], label="Earth")
    ax[0, 1].set_title(f"{dts[1]}")
    ax[0, 1].axis("equal")

    #-----------

    ax[1, 0].plot((E1 - E1[0])/E1[0])
    ax[1, 0].plot((E2 - E2[0])/E2[0])

    #-----------
    
    ax[1, 1].plot(np.arange(len(diff1)), diff1, label="diff1")
    ax[1, 1].plot(np.arange(len(diff2)), diff2, label="diff2")
    ax[1, 1].set_title(f"{dts[3]}")
    plt.legend(loc="best")
    plt.show()

    #-----------

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
    
    # Forward Euler
    #-----------

    filepath = f"data_files/task_5c_fe_dt={dts[0]}.txt"
    data = np.loadtxt(filepath, unpack=True)
    
    tol = 3.5e-3
    idx = np.where( np.abs(data[2]) < tol )
    idxidx = np.where( np.diff(idx[0]) > 10 )
    idx = idx[0][idxidx]
    diff1 = data[1][idx][np.where( (data[1][idx]) > 0 )]

    v = np.sqrt(data[4]**2 + data[5]**2 + data[6]**2)*c.AU*1000/31556926
    r = np.sqrt(data[1]**2 + data[2]**2 + data[3]**2)*c.AU*1000
    K = 1/2*earth_mass*v**2
    V = -G*solar_mass*earth_mass/r
    E = K + V

    # print(np.diff(E)[::100]/E[0])
    print((E[::100] - E[0])/E[0])
    
    ax[0, 0].plot(data[1], data[2], label="Earth")
    ax[0, 0].set_title(f"{dts[0]}")
    ax[0, 0].axis("equal")

    #-----------

    filepath = f"data_files/task_5c_fe_dt={dts[1]}.txt"
    data = np.loadtxt(filepath, unpack=True)

    tol = 3.5e-2
    idx = np.where( np.abs(data[2]) < tol )
    idxidx = np.where( np.diff(idx[0]) > 10 )
    idx = idx[0][idxidx]
    diff2 = data[1][idx][np.where( (data[1][idx]) > 0 )]

    ax[0, 1].plot(data[1], data[2], label="Earth")
    ax[0, 1].set_title(f"{dts[1]}")
    ax[0, 1].axis("equal")

    #-----------

    filepath = f"data_files/task_5c_fe_dt={dts[2]}.txt"
    data = np.loadtxt(filepath, unpack=True)
    ax[1, 0].plot(data[1], data[2], label="Earth")
    ax[1, 0].set_title(f"{dts[2]}")
    ax[1, 0].axis("equal")

    #-----------

    ax[1, 1].plot(np.arange(len(diff1)), diff1, label="diff1")
    ax[1, 1].plot(np.arange(len(diff2)), diff2, label="diff2")

    #-----------

    plt.legend(loc="best")
    plt.show()


if __name__ == "__main__":
    task_5c()
    # data = np.loadtxt("data_files/all_planets.txt", unpack=True)

    # data = np.load("data_files/all_planets.npy")

    # planets = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn",
    #     "Uranus", "Neptune", "Pluto"]
    # for i in range(9):
    #     plt.plot(data[1 + i*6, ::10], data[2 + i*6, ::10], label=planets[i])

    # plt.legend(loc="best")
    # plt.axis("equal")
    # plt.show()
    
    # # plt.plot(data[1], data[2], label="mercury")
    # # plt.plot(data[1], data[2], label="venus")
    # # plt.plot(data[1], data[2], label="earth")
    # # plt.plot(data[1], data[2], label="mars")
    # # plt.plot(data[1], data[2], label="jupiter")
    # # plt.plot(data[1], data[2], label="saturn")
    # # plt.plot(data[1], data[2], label="uranus")
    # # plt.plot(data[1], data[2], label="neptune")
    # # plt.plot(data[1], data[2], label="pluto")

    # # plt.plot(data[1], data[2], label="earth")
    # # plt.plot(data[7], data[8], label="jupiter")
    # # plt.plot(data[13], data[14], label="uranus")

