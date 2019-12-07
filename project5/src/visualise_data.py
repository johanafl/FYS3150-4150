import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import constants as c

def total_energy_and_angular_momentum(data):
    """
    Used for task 5c. Generates the total energy and angular momentum
    for the system at every time step.
    """

    earth_mass = 5.972e24
    solar_mass   = 1.9891e30
    G = 6.67e-11

    rvec = data[1:4].transpose()*c.AU*1000
    vvec = data[4:7].transpose()*c.AU*1000/31556926
    
    r = np.sqrt(data[1]**2 + data[2]**2 + data[3]**2)*c.AU*1000
    v = np.sqrt(data[4]**2 + data[5]**2 + data[6]**2)*c.AU*1000/31556926
    K = 1/2*earth_mass*v**2
    V = -G*solar_mass*earth_mass/r

    L = np.cross(rvec, vvec)*earth_mass
    
    return K + V, L.transpose()


def task_5c():
    dts = ["0.001000", "0.010000"]


    # Velocity Verlet
    #-----------
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
    
    filepath = f"data_files/task_5c_vv_dt={dts[0]}.txt"
    data = np.loadtxt(filepath, unpack=True)
    
    tol = float(dts[0])
    idx = np.where( np.abs(data[2]) < tol )
    diff1 = data[1][idx][np.where( (data[1][idx]) > 0 )]

    E1, L1vv = total_energy_and_angular_momentum(data)

    ax[0, 0].plot(data[1], data[2], label="Earth", color="black")
    ax[0, 0].set_xlabel("Position, [AU]", fontsize=20)
    ax[0, 0].set_ylabel("Position, [AU]", fontsize=20)
    ax[0, 0].set_title(f"dt = {float(dts[0])} yr", fontsize=23)
    ax[0, 0].set_xticks([-1, 0, 1])
    ax[0, 0].set_yticks([-1, 0, 1])
    ax[0, 0].tick_params(labelsize=20)
    ax[0, 0].axis("equal")
    ax[0, 0].grid()
    
    #-----------

    filepath = f"data_files/task_5c_vv_dt={dts[1]}.txt"
    data = np.loadtxt(filepath, unpack=True)
    
    tol = 3e-2
    idx = np.where( np.abs(data[2]) < tol )
    diff2 = data[1][idx][np.where( (data[1][idx]) > 0 )]

    E2, L2vv = total_energy_and_angular_momentum(data)

    ax[0, 1].plot(data[1], data[2], label="Earth", color="gray")
    ax[0, 1].set_xlabel("Position, [AU]", fontsize=20)
    ax[0, 1].set_ylabel("Position, [AU]", fontsize=20)
    ax[0, 1].set_title(f"dt = {float(dts[1])} yr", fontsize=23)
    ax[0, 1].set_xticks([-1, 0, 1])
    ax[0, 1].set_yticks([-1, 0, 1])
    ax[0, 1].tick_params(labelsize=20)
    ax[0, 1].axis("equal")
    ax[0, 1].grid()

    #-----------
    end = 2000
    ax[1, 0].plot(np.abs((E1[0:end] - E1[0])/E1[0]), color="black")
    ax[1, 0].plot(np.abs((E2[0:end] - E2[0])/E2[0]), color="gray")
    ax[1, 0].set_xlabel("Number of time steps", fontsize=20)
    ax[1, 0].set_ylabel(r"Rel. energy error, $[E_0]$", fontsize=20)
    # ax[1, 0].set_title("Total energy", fontsize=23)
    ax[1, 0].set_xticks([0, 1000, 2000])
    ax[1, 0].set_yticks([0, 2.5e-6/2, 2.5e-6])
    ax[1, 0].set_yticklabels([0, r"$1.25 \cdot 10^{-6}$", r"$2.5 \cdot 10^{-6}$"])
    ax[1, 0].tick_params(labelsize=20)
    ax[1, 0].grid()

    #-----------
    
    end = 35
    ax[1, 1].plot(np.abs(diff1[0:end] - diff1[0]), "--o", color="black")
    ax[1, 1].plot(np.abs(diff2[0:end] - diff2[0]), "--o", color="gray")
    ax[1, 1].set_xlabel("Number of orbits", fontsize=20)
    ax[1, 1].set_ylabel("Displacement, [AU]", fontsize=20)
    ax[1, 1].tick_params(labelsize=20)
    ax[1, 1].grid() 
    plt.tight_layout(pad=0.5)
    plt.show()


    #-----------

    
    # Forward Euler
    #-----------
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))

    filepath = f"data_files/task_5c_fe_dt={dts[0]}.txt"
    data = np.loadtxt(filepath, unpack=True)
    
    tol = 3.5e-3
    idx = np.where( np.abs(data[2]) < tol )
    idxidx = np.where( np.diff(idx[0]) > 10 )
    idx = idx[0][idxidx]
    diff1 = data[1][idx][np.where( (data[1][idx]) > 0 )]

    E1, L1fe = total_energy_and_angular_momentum(data)

    ax[0, 0].plot(data[1], data[2], label="Earth", color="black")
    ax[0, 0].set_xlabel("Position, [AU]", fontsize=20)
    ax[0, 0].set_ylabel("Position, [AU]", fontsize=20)
    ax[0, 0].set_title(f"dt = {float(dts[0])} yr", fontsize=23)
    ax[0, 0].set_xticks([-2, 0, 2])
    ax[0, 0].set_yticks([-2, 0, 2])
    ax[0, 0].tick_params(labelsize=20)
    ax[0, 0].axis("equal")
    ax[0, 0].grid()

    #-----------

    filepath = f"data_files/task_5c_fe_dt={dts[1]}.txt"
    data = np.loadtxt(filepath, unpack=True)

    tol = 3.5e-2
    idx = np.where( np.abs(data[2]) < tol )
    idxidx = np.where( np.diff(idx[0]) > 10 )
    idx = idx[0][idxidx]
    diff2 = data[1][idx][np.where( (data[1][idx]) > 0 )]

    E2, L2fe = total_energy_and_angular_momentum(data)

    ax[0, 1].plot(data[1], data[2], label="Earth", color="gray")
    ax[0, 1].set_xlabel("Position, [AU]", fontsize=20)
    ax[0, 1].set_ylabel("Position, [AU]", fontsize=20)
    ax[0, 1].set_title(f"dt = {float(dts[1])} yr", fontsize=23)
    ax[0, 1].set_xticks(np.arange(-7.5, 7.5+5, 5))
    ax[0, 1].set_yticks(np.arange(-8, 4+2, 3))
    ax[0, 1].tick_params(labelsize=20)
    ax[0, 1].axis("equal")
    ax[0, 1].grid()

    #-----------

    end = 10000
    ax[1, 0].plot(np.abs((E1[0:end] - E1[0])/E1[0]), color="black")
    ax[1, 0].plot(np.abs((E2[0:end] - E2[0])/E2[0]), color="gray")
    ax[1, 0].set_xlabel("Number of time steps", fontsize=20)
    ax[1, 0].set_ylabel(r"Rel. energy error, $[E_0]$", fontsize=20)
    # ax[1, 0].set_title("Rel. energy error", fontsize=23)
    ax[1, 0].set_xticks(np.arange(0, 1e4+1, 3000))
    ax[1, 0].tick_params(labelsize=20)
    ax[1, 0].grid()

    #-----------
    
    end = 10
    ax[1, 1].plot(np.abs(diff1[0:end] - diff1[0]), "--o", color="black")
    ax[1, 1].plot(np.abs(diff2[0:end] - diff2[0]), "--o", color="gray")
    ax[1, 1].set_xlabel("Number of orbits", fontsize=20)
    ax[1, 1].set_ylabel("Displacement, [AU]", fontsize=20)
    ax[1, 1].tick_params(labelsize=20)
    ax[1, 1].grid() 

    #-----------

    plt.tight_layout(pad=0.5)
    plt.show()
    end = 10000
    L1vv = np.linalg.norm(L1vv, axis=0)
    L2vv = np.linalg.norm(L2vv, axis=0)
    L1fe = np.linalg.norm(L1fe, axis=0)
    L2fe = np.linalg.norm(L2fe, axis=0)

    #-----------

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
    
    ax.plot(L1vv[0:end]/L1vv[0], label=f"VV, dt={float(dts[0])}", color="black")
    ax.plot(L2vv[0:end]/L2vv[0], label=f"VV, dt={float(dts[1])}", color="black",
        linestyle="dashed")
    ax.plot(L1fe[0:end]/L1fe[0], label=f"FE, dt={float(dts[0])}", color="gray")
    ax.plot(L2fe[0:end]/L2fe[0], label=f"FE, dt={float(dts[1])}", color="gray",
        linestyle="dashed")

    ax.tick_params(labelsize=20)
    ax.grid()
    ax.set_xlabel("Number of time steps", fontsize=20)
    ax.set_ylabel(r"Total angular momentum, $[L_0]$", fontsize=20)
    ax.set_xticks(np.arange(0, 1e4+1, 3000))
    ax.set_yticks(np.arange(1, 2.4+0.4, 0.4))

    ax.legend(fontsize=15, loc="upper left")
    plt.show()


def task_5d():
    data = np.loadtxt("data_files/task_5d.txt", unpack=True)

    plt.plot(data[1], data[2], label="earth")
    plt.legend(loc="best")
    plt.axis("equal")
    plt.show()
    

if __name__ == "__main__":
    task_5c()
    # task_5d()
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

