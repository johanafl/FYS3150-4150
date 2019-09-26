import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class VisualizeData:

    def __init__(self):
        """
        Loads data from text file generated from the accompanying C++ program,
        quantum_dots.cpp.
        
        The file starts with three values which are the number
        of eigenvalues per rho_max, the number of different rho_max values, and
        the number of different grid point values, repectively.
        """

        num_eig, num_rho, num_n = np.loadtxt("eigenvalues.txt", max_rows=1)

        self.num_eig = int(num_eig)      # number of eigenvalues
        self.num_rho = int(num_rho)      # number of rho_max values
        self.num_n   = int(num_n)        # number of grid point values

        self.calc, self.exact, self.error, self.rho_max, self.n = \
            np.loadtxt("eigenvalues.txt", skiprows=2, unpack=True)

    def contour_plot(self):
        """
        Generates contour plot of grid points vs rho_max, and the relative error
        as height values.
        """

        rho_max = np.zeros(self.num_rho)
        error = np.zeros((self.num_n, self.num_rho))
        n = np.zeros(self.num_n)

        for j in range(self.num_n):
            # reading each batch of values per grid size
            
            for i in range(self.num_rho):
                # reading each batch of eigenval errors per rho_max
                start = (408*j) + i*self.num_eig
                stop  = (408*j) + (i + 1)*self.num_eig
                
                idx = np.argmax(self.error[start:stop])

                error[j, i] = np.log10(self.error[start:stop][idx])
                rho_max[i]  = self.rho_max[start:stop][idx]
                n[j] = self.n[start:stop][idx]

        fig, ax = plt.subplots()

        rho_max, n = np.meshgrid(rho_max, n)
        plt.contourf(rho_max, n, error)
        plt.colorbar()

        plt.show()


    def comparison_plot(self):
        """
        Generates a standard plot of the different graphs for comparison. Picks
        the maximum error value per rho_max per grid point value.
        """

        max_rhos  = np.zeros(self.num_rho)
        max_error = np.zeros(self.num_rho)

        fig, ax = plt.subplots()
        
        col = plt.cm.winter(np.linspace(0, 1, self.num_n)) # 9 = max lenght of j
        for j in range(self.num_n):
            # reading each batch of values per grid size

            for i in range(self.num_rho):
                # reading each batch of eigenval errors per rho_max
                start = (408*j) + i*self.num_eig
                stop  = (408*j) + (i + 1)*self.num_eig
                
                idx = np.argmax(self.error[start:stop])

                max_error[i] = self.error[start:stop][idx]
                max_rhos[i]  = self.rho_max[start:stop][idx]
            
        
            ax.plot(max_rhos, max_error, label=100+(j + 1)*10, alpha=0.8, color=col[j])
        
        ax.set_xlabel("exact")
        ax.set_ylabel("error")
        plt.legend()
        plt.tight_layout(pad=2)
        plt.grid()
        plt.show()



if __name__ == "__main__":
    q = VisualizeData()
    q.contour_plot()
    #q.comparison_plot()
    pass