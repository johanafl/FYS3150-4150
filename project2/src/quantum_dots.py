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



def exact_freq_0_05(rho):
    """
    Analytical solution for the eigenvector corresponding to freq = 0.05
    corresponds to n = 3 in article by M. Taut. Eigenvalue (proportional to
    energy) lambda = 0.1750.

    The result is not normalized.

    Parameters
    ----------
    rho : float, numpy.ndarray
        Single value or values of the scaled dimensionless distance, rho.
    """

    l = 0

    return rho**(l + 1)*np.exp(-rho**2/(8*(4*l + 5))) * (1 + rho/(2*(l + 1)) + \
        rho**2/(4*(l + 1)*(4*l + 5)))


def exact_freq_0_25(r):
    """
    Analytical solution for the eigenvector corresponding to freq = 0.25
    corresponds to n = 2 in article by M. Taut. Eigenvalue (proportional to
    energy) lambda = 0.6250.

    The result is not normalized.

    Parameters
    ----------
    rho : float, numpy.ndarray
        Single value or values of the scaled dimensionless distance, rho.
    """

    l = 0

    return r**(l + 1) * np.exp(-r**2/(8*(l + 1))) * (1 + r/(2*(l + 1)))


def approximate_eigenvalues(freq):
    """
    Approximating eigenvalues for the two electron interaction problem.

    See M. Taut. https://journals.aps.org/pra/pdf/10.1103/PhysRevA.48.3561.

    Parameters
    ----------
    freq : float, numpy.ndarray
        Harmonic oscillator frequency.

    Returns
    -------
    : float, numpy.ndarray
        Approximated eigenvalue.

    CURRENTLY NOT IN USE.
    """
    V0     = 3/2*(freq/2)**(2/3)
    freq_e = np.sqrt(3)*freq
    m = 0

    return V0 + freq_e*(m + 1/2)


def approximate_eigenvectors(freq, rho):
    """
    Approximating eigenvectors for the two electron interaction problem.

    See M. Taut. https://journals.aps.org/pra/pdf/10.1103/PhysRevA.48.3561.

    Parameters
    ----------
    freq : float
        Harmonic oscillator frequency.

    rho : numpy.ndarray
        Array with values in the interval [rho_min, rho_max].

    Returns
    -------
    : float, numpy.ndarray
        Approximated eigenvalue.

    CURRENTLY NOT IN USE.
    """

    rho_0  = (2*freq**2)**(-1/3)
    freq_e = np.sqrt(3)*freq

    return (freq_e/np.pi)**(1/4)*np.exp(-(1/2)*freq_e*(rho - rho_0))


def visualize_eigendata_two_electrons_numerical_and_analytical():
    """
    Reads eigenvalue file which contains two columns. The first column consists
    of rho max values, and the second column of the corresponding eigenvalues.

    Reads eigenvector file which contains one column with rho max values and
    then k more columns which contains the k'th element of each eigenvector.

    Plots the two frequencies with analytical results and compares the numerical
    and analytical results. freq = [0.05, 0.25].
    """

    filenames_1 = ["eigenvalue_omega_0.050000.txt", "eigenvalue_omega_0.250000.txt"]
    filenames_2 = ["eigenvector_omega_0.050000.txt", "eigenvector_omega_0.250000.txt"]
    functions = [exact_freq_0_05, exact_freq_0_25]
    exact_eigenvalues = [2*0.1750, 2*0.6250]
    frequencies = [0.05, 0.25]

    for i in range(2):

        _, ax = plt.subplots()

        # eigenvalue calculations
        rho, eigenvalue = np.loadtxt(filenames_1[i], skiprows=1, unpack=True)
        
        error = np.abs(eigenvalue - exact_eigenvalues[i])
        min_error_idx = np.argmin(error)
        min_error = error[min_error_idx]
        
        # recheck == True if best rho max is at the start or end of interval
        recheck = (min_error_idx == (len(rho) - 1)) or (min_error_idx == 0)
        
        print("eigenvalue data")
        print("===============")
        print("filename: ", filenames_1[i])
        print("exact eigenvalue: ", exact_eigenvalues[i])
        print("best eigenvalue: ", eigenvalue[min_error_idx])
        print("best error: ", min_error)
        print("best rho max: ", rho[min_error_idx])
        print("index of rho max: ", min_error_idx)
        print("numbero rho maxo: ", len(rho))
        print("recheck: ", recheck, "\n")


        # eigenvector calculations

        eigenvectors     = np.loadtxt(filenames_2[i], skiprows=1)
        num_rho_max      = np.shape(eigenvectors)[0]        # number of rho max values
        num_eig_elements = np.shape(eigenvectors)[1] - 1    # number of elements in eigenvector

        rhos, eigenvectors = eigenvectors[:, 0], eigenvectors[:, 1:]
        eigenvectors[min_error_idx] /= np.sum(eigenvectors[min_error_idx]) # normalizing

        rho = np.linspace(0, rhos[min_error_idx], num_eig_elements)
        exact_eigenvector  = functions[i](rho)
        exact_eigenvector /= np.sum(exact_eigenvector)
        
        ax.plot(rho, exact_eigenvector, label="exact")
        ax.plot(rho, eigenvectors[min_error_idx], label="computed")
        ax.set_xlabel(r"$\rho$", fontsize=40)
        ax.set_ylabel(r"$u(\rho)$", fontsize=40)
        ax.set_title(r"$\omega_r = $" + f"{frequencies[i]}", fontsize=40)
        
        ax.tick_params(labelsize=30)
        ax.grid()
        plt.legend(fontsize=30)
        plt.show()


def visualize_eigendata_two_electrons_numerical():
    """
    Reads eigenvalue file which contains two columns. The first column consists
    of rho max values, and the second column of the corresponding eigenvalues.

    Reads eigenvector file which contains one column with rho max values and
    then k more columns which contains the k'th element of each eigenvector.

    Presents the radial part of the wave function for all the different
    frequencies.
    """

    filenames_1 = ["eigenvalue_omega_0.010000.txt", "eigenvalue_omega_0.050000.txt", 
        "eigenvalue_omega_0.250000.txt", "eigenvalue_omega_0.500000.txt",
        "eigenvalue_omega_1.000000.txt", "eigenvalue_omega_5.000000.txt"]
    filenames_2 = ["eigenvector_omega_0.010000.txt", "eigenvector_omega_0.050000.txt",
        "eigenvector_omega_0.250000.txt", "eigenvector_omega_0.500000.txt",
        "eigenvector_omega_1.000000.txt", "eigenvector_omega_5.000000.txt"]

    frequencies = [0.01, 0.05, 0.25, 0.5, 1, 5]


    _, ax = plt.subplots()

    for i in range(6):

        # eigenvalue calculations
        rho, eigenvalue = np.loadtxt(filenames_1[i], skiprows=1, unpack=True)
        
        best_rho_idx = np.argmin(np.abs(eigenvalue[:-1] - eigenvalue[1:]))

        best_rho = rho[best_rho_idx]

        print("eigenvalue data")
        print("===============")
        print("filename: ", filenames_1[i])
        print("best rho: ", best_rho)
        print("best rho idx: ", best_rho_idx)
        print("number of rhos: ", len(rho))
        


        # eigenvector calculations

        eigenvectors     = np.loadtxt(filenames_2[i], skiprows=1)
        num_rho_max      = np.shape(eigenvectors)[0]        # number of rho max values
        num_eig_elements = np.shape(eigenvectors)[1] - 1    # number of elements in eigenvector

        rhos, eigenvectors = eigenvectors[:, 0], eigenvectors[:, 1:]
        eigenvectors[best_rho_idx] /= np.sum(eigenvectors[best_rho_idx]) # normalizing

        rho = np.linspace(0, rhos[best_rho_idx], num_eig_elements)

    
        ax.plot(rho, eigenvectors[best_rho_idx], label=f"$\omega_r: ${frequencies[i]}")
        ax.set_xlabel(r"$\rho$", fontsize=40)
        ax.set_ylabel(r"$u(\rho)$", fontsize=40)
    
    ax.grid()
    ax.tick_params(labelsize=30)
    plt.legend(fontsize=30)
    plt.show()



if __name__ == "__main__":
    q = VisualizeData("eigenvalues.txt")
    q.contour_plot(selection="max")
    q.comparison_plot()

    visualize_eigendata_two_electrons_numerical_and_analytical()
    visualize_eigendata_two_electrons_numerical()
    pass