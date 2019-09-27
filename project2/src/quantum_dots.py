import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class VisualizeData:

    def __init__(self, filename="eigenvalues_tmp.txt"):
        """
        Loads data from text file generated from the accompanying C++ program,
        quantum_dots.cpp.
        
        The file starts with three values which are the number
        of eigenvalues per rho_max, the number of different rho_max values, and
        the number of different grid point values, repectively.
        """

        num_eig, num_rho, num_n = np.loadtxt(filename, max_rows=1)

        self.num_eig = int(num_eig)      # number of eigenvalues
        self.num_rho = int(num_rho)      # number of rho_max values
        self.num_n   = int(num_n)        # number of grid point values

        self.calc, self.exact, self.error, self.rho_max, self.n = \
            np.loadtxt(filename, skiprows=2, unpack=True)


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
                start = self.num_eig*(self.num_rho + 1)*j + i*self.num_eig
                stop  = self.num_eig*(self.num_rho + 1)*j + (i + 1)*self.num_eig
                
                idx = np.argmax(self.error[start:stop])

                error[j, i] = np.log10(self.error[start:stop][idx])
                rho_max[i]  = self.rho_max[start:stop][idx]
                n[j] = self.n[start:stop][idx]

        fig, ax = plt.subplots()

        rho_max, n = np.meshgrid(rho_max, n)
        plt.contourf(rho_max, n, error)
        plt.colorbar()

        plt.show()


    def comparison_plot(self, ax=None):
        """
        Generates a standard plot of the different graphs for comparison. Picks
        the maximum error value per rho_max per grid point value.
        """

        max_rhos  = np.zeros(self.num_rho)
        max_error = np.zeros(self.num_rho)

        if ax is None:
            _, ax = plt.subplots()
        
        col = plt.cm.winter(np.linspace(0, 1, self.num_n))
        for j in range(self.num_n):
            # reading each batch of values per grid size

            for i in range(self.num_rho):
                # reading each batch of eigenval errors per rho_max
                start = self.num_eig*(self.num_rho + 1)*j + i*self.num_eig
                stop  = self.num_eig*(self.num_rho + 1)*j + (i + 1)*self.num_eig
                
                idx = np.argmax(self.error[start:stop])

                max_error[i] = self.error[start:stop][idx]
                max_rhos[i]  = self.rho_max[start:stop][idx]
            
            if ax is None:
                ax.plot(max_rhos, max_error, label=100+(j + 1)*10, alpha=0.8, color=col[j])

            else:
                ax.plot(max_rhos, max_error)

        if ax is None:
            ax.set_xlabel("exact")
            ax.set_ylabel("error")
            plt.legend()
            plt.tight_layout(pad=2)
            plt.grid()
            plt.show()

    def visualize_frequencies(self):

        fig, ax = plt.subplots()
        
        filenames = ["eigenvalues_w_0.010000.txt", "eigenvalues_w_0.500000.txt", 
                    "eigenvalues_w_1.000000.txt", "eigenvalues_w_5.000000.txt"]

        for filename in filenames:
            
            with open(filename, "r"):
                self.__init__(filename)
                self.comparison_plot(ax)
        
        ax.set_xlabel("exact")
        ax.set_ylabel("error")
        plt.legend()
        plt.tight_layout(pad=2)
        plt.grid()
        plt.show()



if __name__ == "__main__":

    filenames = ["eigenvalues_w_0.010000.txt", "eigenvalues_w_0.500000.txt", 
                "eigenvalues_w_1.000000.txt", "eigenvalues_w_5.000000.txt"]
    
    # filenames = ["eigenvalues_tmp.txt"]

    for filename in filenames:
        q = VisualizeData(filename)
        q.contour_plot()
    # q = VisualizeData()
    # q.visualize_frequencies()
    pass