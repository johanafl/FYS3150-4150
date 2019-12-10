import time
import os
import os.path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# All masses in kg.
solar_mass     = 1.9891e30
mercury_mass = 3.285e23
venus_mass   = 4.867e24
earth_mass   = 5.972e24
mars_mass    = 6.39e23
jupiter_mass = 1.898e27
saturn_mass  = 5.683e26
uranus_mass  = 8.681e25
neptune_mass = 1.024e26
pluto_mass   = 1.309e22
yr = 31556926   # Seconds in a year.
AU = 149597871*1000 # Meters in one AU.
G = 6.67e-11                # [m^3 kg^-1 s^-2]


def total_energy_and_angular_momentum(data):
    """
    Used for task 5c. Generates the total energy and angular momentum
    for the system at every time step.
    """


    rvec = data[1:4].transpose()*AU
    vvec = data[4:7].transpose()*AU/yr
    
    r = np.sqrt(data[1]**2 + data[2]**2 + data[3]**2)*AU
    v = np.sqrt(data[4]**2 + data[5]**2 + data[6]**2)*AU/yr
    K = 1/2*earth_mass*v**2
    V = -G*solar_mass*earth_mass/r

    L = np.cross(rvec, vvec)*earth_mass
    
    return K + V, L.transpose()

def total_energy_earth_and_jupiter_5e(data):
    
    rvec_earth   = data[1:4]*AU # Sun to Earth.
    rvec_jupiter = data[7:10]*AU # Sun to Jupiter.
    rvec_earth_jupiter = np.abs(rvec_earth - rvec_jupiter)  # Earth to Jupiter.

    vvec_earth = data[4:7]*AU/yr
    vvec_jupiter = data[10:13]*AU/yr

    r_earth = np.linalg.norm(rvec_earth, axis=0)
    r_jupiter = np.linalg.norm(rvec_jupiter, axis=0)
    r_earth_jupiter = np.linalg.norm(rvec_earth_jupiter, axis=0)

    v_earth = np.linalg.norm(vvec_earth, axis=0)
    v_jupiter = np.linalg.norm(vvec_jupiter, axis=0)
    
    K_earth = 1/2*earth_mass*v_earth**2
    K_jupiter = 1/2*jupiter_mass*v_jupiter**2
    
    V_earth = -G*solar_mass*earth_mass/r_earth
    V_jupiter = -G*solar_mass*jupiter_mass/r_jupiter
    V_earth_jupiter = -G*earth_mass*jupiter_mass/r_earth_jupiter

    
    return K_earth + K_jupiter + V_earth + V_jupiter + V_earth_jupiter

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


def task_5d_escape_velocity():
    
    energies = []
    velocities = []
    directory = "data_files/"

    for data_file in os.listdir(directory):
        # Loops over all files in directory.
        
        filename = os.fsdecode(data_file)
        if filename.startswith("task_5d"):
            
            data = np.loadtxt(f"data_files/{filename}", unpack=True)
            E, _ = total_energy_and_angular_momentum(data)
            energies.append(E)
            
            velocities.append(float(filename[8:-4]))
    
    velocities = np.array(sorted(velocities))
    energies = np.array(sorted(energies))

    idx = np.argmin( np.abs(energies) )
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    fig.text(x=0.51, y=0.61, s=f"({velocities[idx]:.1f}, 0)", fontsize=20, color="black")
    
    ax.plot(velocities, energies, color="black")
    ax.plot(velocities[idx], energies[idx], "ko")
    ax.set_xlabel("Initial velocity, [AU/yr]", fontsize=20)
    ax.set_ylabel("Initial total energy, [J]", fontsize=20)
    ax.set_xticks(np.arange(6.5, 10+1, 1))
    ax.set_yticks(np.arange(-2.5e33, 1.5e33+1e33, 1e33))
    ax.tick_params(labelsize=20)
    ax.grid()
    
    plt.show()


def task_5d_beta():

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(11, 9))
    ax = ax.reshape(-1)

    filenames = []
    
    directory = "data_files/"
    
    for data_file in os.listdir(directory):
        # Loops over all files in directory.
        filename = os.fsdecode(data_file)
        
        if filename.startswith("varying_beta"):
            filenames.append(filename)

    filenames = sorted(filenames)
    for i in range(4):
        data = np.loadtxt(f"data_files/{filenames[i]}", unpack=True)
        ax[i].plot(data[1], data[2], color="black")
        ax[i].tick_params(labelsize=20)
        ax[i].set_title(r"$\beta$ = " + f"{float(filenames[i][13:-4])}", fontsize=20)
        ax[i].grid()
        ax[i].axis("equal")



    fig.text(x=0.02, y=0.4, s="Position, [AU]", fontsize=20, rotation="vertical")
    fig.text(x=0.42, y=0.03, s="Position, [AU]", fontsize=20)
    # ax.legend()
    # ax.set_xticks(np.arange(6.5, 10+1, 1))
    # ax.set_yticks(np.arange(-2.5e33, 1.5e33+1e33, 1e33))
    
    plt.show()


def task_5e():
    fig1, ax1 = plt.subplots(nrows=2, ncols=2, figsize=(11, 9))
    fig1.text(x=0.02, y=0.4, s="Position, [AU]", fontsize=20, rotation="vertical")
    fig1.text(x=0.42, y=0.03, s="Position, [AU]", fontsize=20)
    fig2, ax2 = plt.subplots(nrows=2, ncols=2, figsize=(11, 9))
    fig2.text(x=0.02, y=0.31, s="Distance from the Sun, [AU]", fontsize=20, rotation="vertical")
    fig2.text(x=0.46, y=0.03, s="Time, [yr]", fontsize=20)
    fig3, ax3 = plt.subplots(nrows=2, ncols=2, figsize=(11, 9))
    fig3.text(x=0.02, y=0.38, s=r"Total energy, $[E_0]$", fontsize=20, rotation="vertical")
    fig3.text(x=0.46, y=0.03, s="Time, [yr]", fontsize=20)
    ax1 = ax1.reshape(-1)
    ax2 = ax2.reshape(-1)
    ax3 = ax3.reshape(-1)

    dt = 1e-5   # Time step length, used for scaling the x axis.
    filenames = []
    masses = []
    
    directory = "data_files/"
    
    for data_file in os.listdir(directory):
        # Loops over all files in directory.
        filename = os.fsdecode(data_file)
        
        if filename.startswith("jupiter_mass") and filename.endswith(".txt"):
            filenames.append(filename)
            masses.append(float(filename[13:-4]))

    filenames = [x for _, x in sorted(zip(masses, filenames))]


    
    for i in range(4):
        try:
            data = np.load("data_files/" + filenames[i][:-4] + ".npy")

        except FileNotFoundError:
            convert_to_npy("data_files/" + filenames[i][:-4])
            data = np.load("data_files/" + filenames[i][:-4] + ".npy")


        ax1[i].plot(data[1], data[2], label="Earth")
        ax1[i].plot(data[7], data[8], label="Jupiter")
        ax1[i].tick_params(labelsize=20)
        ax1[i].set_title(r"$M_{Jupiter}$ = " + f"{float(filenames[i][13:-4]):.3e}", fontsize=20)
        ax1[i].grid()
        ax1[i].axis("equal")

        earth_dist = np.linalg.norm(data[1:4].transpose(), axis=1)
        ax2[i].plot(np.arange(1, len(earth_dist)+1, 1)*dt, earth_dist, color="black")
        ax2[i].tick_params(labelsize=20)
        ax2[i].set_title(r"$M_{Jupiter}$ = " + f"{float(filenames[i][13:-4]):.3e}", fontsize=20)
        ax2[i].grid()

        E = total_energy_earth_and_jupiter_5e(data)
        ax3[i].plot(np.arange(1, len(E)+1, 1)*dt, E/E[0], color="black")
        ax3[i].tick_params(labelsize=15)
        ax3[i].set_title(r"$M_{Jupiter}$ = " + f"{float(filenames[i][13:-4]):.3e}", fontsize=20)
        # ax3[i].yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        ax3[i].grid()

    plt.show()


def task_5f():

    data = np.loadtxt(f"data_files/sun_earth_jupiter.txt", unpack=True)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    ax.plot(data[7], data[8], label="Earth")
    ax.plot(data[13], data[14], label="Jupiter")
    ax.plot(data[1], data[2], label="Sun")
    ax.set_xlabel("Position, [AU]", fontsize=20)
    ax.set_ylabel("Position, [AU]", fontsize=20)
    ax.tick_params(labelsize=20)
    ax.grid()
    ax.axis("equal")
    ax.legend(fontsize=20)
    
    plt.show()


def task_5f_all_planets():
    """
    NB: REMEMBER TO DELETE THE .npy FILE IF THE .txt DATA FILE IS UPDATED.
    ELSE YOU WILL JUST VIEW THE SAME DATA OVER AND OVER.
    """
    filename = "data_files/all_planets"
    
    try:
        data = np.load(filename + ".npy")

    except FileNotFoundError:
        convert_to_npy(filename)
        data = np.load(filename + ".npy")
    
    fig, ax = plt.subplots(figsize=(11, 9))

    planets = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn",
        "Uranus", "Neptune", "Pluto"]
    for i in range(10):
        ax.plot(data[1 + i*6, ::10], data[2 + i*6, ::10], label=planets[i], linewidth=0.8)

    
    ax.set_xlabel("Position, [AU]", fontsize=20)
    ax.set_ylabel("Position, [AU]", fontsize=20)
    ax.tick_params(labelsize=20)
    ax.grid()
    ax.axis("equal")
    ax.legend(fontsize=15, loc="upper left")
    plt.show()


def convert_to_npy(filename):
    """
    For converting .txt files to .npy.

    Parameters
    ----------
    filename : str
        Filepath and name without file extension. If .txt file extension
        is present, it will be removed.
    """

    if filename[-4:] == ".txt":
        filename = filename[:-4]    # Removing extension.

    print(f"Converting {filename}.txt to Numpy binary...")
    t1 = time.time()

    data = np.loadtxt(filename + ".txt", unpack=True)
    np.save(filename + ".npy", data)

    print(f"Numpy binary saved to {filename}.npy in {time.time() - t1:.4f} seconds.")

if __name__ == "__main__":
    # task_5c()
    # task_5d_escape_velocity()
    # task_5d_beta()
    task_5e()
    # task_5f()
    # task_5f_all_planets()
    pass

