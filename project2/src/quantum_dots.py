import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import glob, os

class VisualizeData:

    def __init__(self, filename="eigenvalues_tmp.txt"):
        """
        Loads data from text file generated from the accompanying C++ program,
        quantum_dots.cpp.
        
        The file starts with three values which are the number
        of eigenvalues per rho_max, the number of different rho_max values, and
        the number of different grid point values, repectively.
        """

        try:
            num_eig, num_rho, num_grid = np.loadtxt(filename, max_rows=1)
            self.calc, self.exact, self.error, self.rho_max, self.grid = \
                np.loadtxt(filename, skiprows=2, unpack=True)

        except:
            msg = (f"Could not read file with filename {filename}.\n"
                "Perhaps this is not a file with eigenvalues/errors for a "
                "quantum mechanical system?")
            raise TypeError(msg)

        self.num_eig  = int(num_eig)    # number of eigenvalues
        self.num_rho  = int(num_rho)    # number of rho_max values
        self.num_grid = int(num_grid)   # number of grid point values

    def contour_plot(self, selection="max"):
        """
        Generates contour plot of grid points vs rho_max, and the relative error
        as height values.

        Parameters
        ----------

        selection : str, int
            For choosing which selection of eigenvalues to visualize. 'min' for
            displaying the minimum error for each rho and grid value, 'max' for
            the maximum error. 0, 1, ... , 7 for displaying error of a specific
            eigenvalue/energy.

        Raises
        ------
        ValueError
            If selection input is an invalid input.
        """

        allowed_selections = [0, 1, 2, 3, 4, 5, 6, 7, "min", "max"]

        if selection not in allowed_selections:
            msg = (f"Selection {selection} not an allowed input. Use one of "
            f"{allowed_selections}.")
            raise ValueError(msg)

        # arrays for storing extracted values
        rho_max = np.zeros(self.num_rho)
        error   = np.zeros((self.num_grid, self.num_rho))
        grid    = np.zeros(self.num_grid)

        for j in range(self.num_grid):
            # reading each batch of values per grid size
            for i in range(self.num_rho):
                # reading each batch of eigenval errors per rho_max

                # generating indices for slicing.
                start = self.num_eig*(self.num_rho + 1)*j + i*self.num_eig
                stop  = self.num_eig*(self.num_rho + 1)*j + (i + 1)*self.num_eig
                
                if selection == "max":
                    # finding the index of the max error per rho_max
                    idx = np.argmax(self.error[start:stop])

                elif selection == "min":
                    # finding the index of the min error per rho_max
                    idx = np.argmin(self.error[start:stop])

                else:
                    # choosing a specific eigenvalue/energy
                    idx = selection

                error[j, i] = np.log10(self.error[start:stop][idx])
                rho_max[i]  = self.rho_max[start:stop][idx]
                grid[j]     = self.grid[start:stop][idx]

        rho_max, grid = np.meshgrid(rho_max, grid)
        
        fig, ax = plt.subplots(figsize=(10, 8))
        plt.contourf(rho_max, grid, error)
        
        cbar = plt.colorbar()
        cbar.set_label(r"$log_{10}$ error", fontsize=40)
        cbar.ax.tick_params(labelsize=30) 

        ax.set_xlabel(r"$\rho_{max}$", fontsize=40)
        ax.set_ylabel(r"grid, n", fontsize=40)
        ax.tick_params(labelsize=30)

        plt.show()


    def comparison_plot(self):
        """
        Generates a standard plot of three different graphs for comparison.
        Picks the maximum error value per rho_max per grid point value.
        """

        max_rhos  = np.zeros(self.num_rho)
        max_error = np.zeros(self.num_rho)

        _, ax = plt.subplots(figsize=(10, 8))

        selected_values = [0, int(self.num_grid/2), self.num_grid-1]
        colors = ["silver", "gray", "black"]
        
        # for j in range(self.num_grid):
        for k, j in enumerate(selected_values):
            # reading each batch of values per grid size

            for i in range(self.num_rho):
                # reading each batch of eigenval errors per rho_max

                # generating indices for slicing.
                start = self.num_eig*(self.num_rho + 1)*j + i*self.num_eig
                stop  = self.num_eig*(self.num_rho + 1)*j + (i + 1)*self.num_eig
                
                # finding the index of the max error per rho_max
                idx = np.argmax(self.error[start:stop])

                # extracts the max error and corresponding rho_max
                max_error[i] = self.error[start:stop][idx]
                max_rhos[i]  = self.rho_max[start:stop][idx]

            
            ax.plot(max_rhos, max_error, color=colors[k],
                label=f"n: {self.grid[start:stop][0]:.0f}", linewidth=3)



        ax.tick_params(labelsize=30)
        ax.set_xlabel(r"$\rho_{max}$", fontsize=40)
        ax.set_ylabel("error", fontsize=40)
        
        plt.legend(fontsize=35)
        plt.grid()
        plt.show()



if __name__ == "__main__":
    q = VisualizeData("eigenvalues.txt")
    # q.contour_plot(selection="min")
    q.comparison_plot()
    pass