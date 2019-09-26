import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def visualize_eigenvalues():
    """
    Visualizing calculated, exact and error.
    """

    num_eig, num_rho = np.loadtxt("eigenvalues.txt", max_rows=1)

    num_eig = int(num_eig)      # number of eigenvalues
    num_rho = int(num_rho)      # number of rho_max values

    max_error = np.zeros(num_rho)
    max_rhos  = np.zeros(num_rho)
    
    calc, exact, error, rho_max, n = \
        np.loadtxt("eigenvalues.txt", skiprows=2, unpack=True)


    x = np.zeros(num_rho)
    y = np.zeros(9)
    Z = np.zeros((9, num_rho))

    for j in range(9):
        for i in range(num_rho):
            # reading each batch of eigenval errors per rho_max
            start = (408*j) + i*num_eig
            stop  = (408*j) + (i + 1)*num_eig
            
            idx = np.argmax(error[start:stop])

            Z[j, i] = np.log10(error[start:stop][idx])
            x[i]    = rho_max[start:stop][idx]
            y[j]    = n[start:stop][idx]

    fig, ax = plt.subplots()
    
    # plt.set_cmap("gist_gray")

    X, Y = np.meshgrid(x, y)
    plt.contourf(X, Y, Z)
    plt.colorbar()

    plt.show()

    
    col = plt.cm.winter(np.linspace(0, 1, 9)) # 9 = max lenght of j
    for j in range(9):
        # reading each batch of values per grid size

        for i in range(num_rho):
            # reading each batch of eigenval errors per rho_max
            start = (408*j) + i*num_eig
            stop  = (408*j) + (i + 1)*num_eig
            
            idx = np.argmax(error[start:stop])

            max_error[i] = error[start:stop][idx]
            max_rhos[i]  = rho_max[start:stop][idx]
        
    
        ax.plot(max_rhos, max_error, label=100+(j + 1)*10, alpha=0.8, color=col[j])
    
    ax.set_xlabel("exact")
    ax.set_ylabel("error")
    plt.legend()
    plt.tight_layout(pad=2)
    plt.grid()
    plt.show()


if __name__ == "__main__":
    visualize_eigenvalues()
    pass