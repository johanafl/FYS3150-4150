import matplotlib.pyplot as plt
import numpy as np
import sys, os
from progress.bar import IncrementalBar

def plot_func(axis, n, args=False, exact=False, filename=False, method=True):

    try:
        if not args:
            values = np.loadtxt(filename.format(n), skiprows=1)
    except ValueError:
        values = args
    
    if method:

        x     = np.linspace(0, 1, n+2)
        u     = values[:, 0]
        v     = values[:, 1]
        error = values[:, 2]

        axis[0].plot(x, v, "-", label=method.format(n))
        axis[1].loglog(n, np.max(error), "o", label="exact, n={:d}".format(n))

        if exact:
            axis[0].plot(x, u, label="u(x)")
    
    else:
        epsilon_general = args
        epsilon_special = filename

        axis[0].plot(n, epsilon_general, label=r"Max(\epsilon(n)),\
                                                 general solution.")
        axis[0].set_xlabel("n")
        axis[0].set_ylabel("Relative error")
        axis[1].plot(n, epsilon_special, label=r"Max(\epsilon(n)),\
                                                 special solution.")
        axis[1].set_xlabel("n")
        axis[1].set_ylabel("Relative error")

if __name__ == "__main__":
    num_of_args = len(sys.argv)
    fig, ax  = plt.subplots(1, 2)

    arg = np.array(["Methheds", "Epstasy"])
    run = arg[1]
    
    # if num_of_args >= 2:
    #     num_n    = num_of_args - 1
    #     exact    = False

    #     for i in range(1, num_n + 1):
    #         if i == num_n:
    #             exact = True

    #         filename = str(sys.argv[i])
    #         values   = np.loadtxt(filename, skiprows=1)
    #         n        = len(values[:, 0]) - 2

    #         plot_func(ax, n, args=values, exact=exact)
        
    #     ax[0].legend()
    #     plt.show()

    if run == arg[0]:
        filenames = np.array(["Poisson_values_n_{:d}.txt",
                              "Poisson_values_spes_n_{:d}.txt",
                              "Poisson_values_LU_n_{:d}.txt"])
        methods   = np.array(["General Thomas n={:d}", "Special Thomas n={:d}",
                              "LU armadillo n={:d}"])

        num_n = 3
        exact = False

        for j in range(3):
            exact=False
            for i in range(1, num_n + 1):
                if i == num_n:
                    exact = True
                plot_func(ax, 10**i, exact=exact, filename=filenames[j],
                            method=methods[j])
    
            ax[0].legend()
            plt.savefig(methods[j].format(1000))
    
    else:
        filenames = np.array(["Poisson_values_n_{:d}.txt",
                              "Poisson_values_spes_n_{:d}.txt",
                              "Poisson_values_LU_n_{:d}.txt"])
        methods   = np.array(["General Thomas n={:d}", "Special Thomas n={:d}",
                              "LU armadillo n={:d}"])
        
        num_n      = 71
        n_array    = np.linspace(1e1, 1e7, 61)
        error_gen  = np.zeros(num_n-10)
        error_spes = np.zeros(num_n-10)

        # bar = IncrementalBar("Processing", max=61)

        # errors = open("relative_errors.txt", "w")
        # errors.write("Thomas              |   Thomas special\n")

        # for i in range(10, num_n):
        #     # str = "./proj1_test.out {:d}".format(int(10**(i/10.0)))
        #     # os.system(str)
        #     values_gen    = np.loadtxt("Poisson_values_n_{:d}.txt"
        #                             .format(int(10**(i/10.0))), skiprows=1)
        #     values_spes   = np.loadtxt("Poisson_values_spes_n_{:d}.txt"
        #                             .format(int(10**(i/10.0))), skiprows=1)
        #     error_gen[i-10]  = np.max(values_gen[:, 2])
        #     error_spes[i-10] = np.max(values_spes[:, 2])
        #     errors.write(str(error_gen[i-10]) + "   " + str(error_spes[i-10])
        #                  + "\n")

        #     bar.next()
        # bar.finish()
        # errors.close()
        
        plot_func(ax, n_array, args=error_gen, filename=error_spes,
                     method=False)
    
        ax[0].legend(loc="best")
        ax[1].legend(loc="best")
        plt.savefig(r"\epsilon(n) (Maximum relative error per n)")