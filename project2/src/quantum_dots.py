import numpy as np
import matplotlib.pyplot as plt

class VisualizeData:
    """
    For visualizing data from task 2d.
    """

    def __init__(self, filename="eigenvalues_tmp.txt"):
        """
        Loads data from text file generated from the accompanying C++ program,
        quantum_dots.cpp.
        
        The file starts with three values which are the number
        of eigenvalues per rho_max, the number of different rho_max values, and
        the number of different grid point values, repectively.

        Parameters
        ----------
        filename : str
            Filename of file to be read. Defaults to 'eigenvalues_tmp.txt'.

        Raises
        ------
        ValueError
            Raised if input file is of incorrect formatting or not existing.
        """

        try:
            num_eig, num_rho, num_grid = np.loadtxt(filename, max_rows=1)
            self.calc, self.exact, self.error, self.rho_max, self.grid = \
                np.loadtxt(filename, skiprows=2, unpack=True)

        except:
            msg = (f"Could not read file with filename {filename}.\n"
                "Perhaps this is not a file with eigenvalues/errors for a "
                "quantum mechanical system?")
            raise ValueError(msg)

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
                # idx = np.argmax(self.error[start:stop])
                idx = 0

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





# Exact radial wavefunction for two electrons. (NB! this is u(r), so not 
# actually the wavefunction.)
# Need r_vector (same as rho). Must be n long.
# We calculate for l = 0, so set l = 0.
# Need an n to loop over/set resolution.

def exact_omega_0_05(r,l=0):
    """
    omega_r = 0.25 (corresponds to n = 2 in article by M. Taut.) -> eigenvalue (proportional to energy) lambda = 0.6250
    """
    return r**(l+1) * np.exp(-r**2/(8*(4*l+5))) * (1 + r/(2*(l+1)) + r**2/(4*(l+1)*(4*l+5))) # NB!!!!! Not normalized!
#     # u[i] = pow(r,l+1) * exp(-r*r/(8*(4*l + 5))) * (1 + r/(2*(l + 1)) + r*r/(4*(l + 1)*(4*l + 5)));

def exact_omega_0_25(r,l=0):
    """
    omega_r = 0.05 (corresponds to n = 3 in article by M. Taut.) -> eigenvalue (proportional to energy) lambda = 0.1750
    """
    return r**(l+1) * np.exp(-r**2/(8*(l+1))) * (1 + r/(2*(l+1))) # NB!!!!! Not normalized!
    # u[i] = pow(r,l+1) * exp(-r*r/(8*(l + 1))) * (1 + r/(2*(l + 1)));


def visualize_eigenvalue_data_two_electrons():
    filename = "eigenvalue_omega_0.250000.txt"
    exact_eigenvalue = 2*0.6250
    rho, eigenvalue = np.loadtxt(filename, skiprows=1, unpack=True)
    min_error = np.min(np.abs(eigenvalue-exact_eigenvalue))
    min_error_idx = np.argmin(np.abs(eigenvalue-exact_eigenvalue))
    # print(min_error_idx)
    # print(eigenvalue)
    print(eigenvalue[min_error_idx])
    print(min_error)
    print(rho)
    print(len(rho))
    return min_error_idx

def visualize_eigenvector_data_two_electrons(idx):
    filename = "eigenvector_omega_0.250000.txt"
    eigenvectors = np.loadtxt(filename, skiprows=1)
    num_rho_max = np.shape(eigenvectors)[0]
    num_eig_elements = np.shape(eigenvectors)[1] - 1

    rhos, eigenvectors = eigenvectors[:, 0], eigenvectors[:, 1:]

    rho = np.linspace(0, rhos[idx], num_eig_elements)
    plt.plot(rho, exact_omega_0_05(rho))
    # u_of_rho = exact_omega_0_25(rho)
    # plt.plot(rho, u_of_rho/np.trapz(u_of_rho),label="exact")
    # plt.plot(rho, eigenvectors[idx],label="computed")
    plt.xlabel(r"$\rho$")
    plt.ylabel(r"eigenvector")
    plt.legend()
    plt.show()

    # for i in range(num_rho_max):
    #     rho_interval = np.linspace(0, rho[i], num_eig_elements)

    #     plt.plot(rho_interval, eigenvectors[i])
    #     plt.xlabel(r"$\rho$")
    #     plt.ylabel(r"eigenvector")
    #     plt.show()


if __name__ == "__main__":
    # q = VisualizeData("eigenvalues.txt")
    # q.contour_plot(selection="max")
    # q.comparison_plot()

    idx = visualize_eigenvalue_data_two_electrons()
    visualize_eigenvector_data_two_electrons(idx)
    pass