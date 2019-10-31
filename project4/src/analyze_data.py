import numpy as np
import matplotlib.pyplot as plt

def analayzan(ax, prop, test=True, mean_E=False, mean_M=False, spec_heat=False, susceptibility=False):
    """
    Visualizes the mean energy, E, mean absolute magnetization, |M|,
    the specific heat, Cv, and the susceptibility, X.

    Parameters
    ----------
    ax : 

    file : 
    """
    E  = file[0]
    M  = np.abs(file[1])
    MC = file[some_index]

    if test:
        if mean_E:
            ax.set_title("Energy of crystal with zomezhings")
            ax.plot(MC, E)
            ax.set_ylabel("E")
        
        if mean_M:
            ax.set_title("Magnetization of crystal with zomezhings")
            ax.plot(MC, M)
            ax.set_ylabel("|M|")
        
        if sepc_heat:
            ax.set_title("Specific heat for ze crystal with zomezhings")
            ax.plot(MC, Cv)
            ax.set_ylabel("$C_V$")
        
        if susceptibility:
            ax.set_title("SUSCEPTUBALATAAAA of crystal with zomezhings")
            ax.plot(MC, X)
            ax.set_ylabel(r"$\Chi$")
    
    else:
        if mean_E:
            ax.set_title("Energy of crystal with zomezhings")
            ax.plot(MC, E, label=r"$T = $%d" %temp)
            ax.set_ylabel("E")
        
        if mean_M:
            ax.set_title("Magnetization of crystal with zomezhings")
            ax.plot(MC, M, label=r"$T = $%d" %temp)
            ax.set_ylabel("|M|")
        
        if sepc_heat:
            ax.set_title("Specific heat for ze crystal with zomezhings")
            ax.plot(MC, Cv, label=r"$T = $%d" %temp)
            ax.set_ylabel("$C_V$")
        
        if susceptibility:
            ax.set_title("SUSCEPTUBALATaaaaa of crystal with zomezhings")
            ax.plot(MC, X, label=r"$T = $%d" %temp)
            ax.set_ylabel(r"$\Chi$")


if __name == "__main__":
    bools   = np.array([True, False])
    fig, ax = plt.subplots()

    test = bools[0]

    if test:
        boolians   = np.array([False, False, False, True, False, False, False])
        properties = np.array(["E", "M", "Cv", "X"])

        for i, element in enumerate(properties):
            file_ = np.loadtxt("{: s}.txt".format(element))

            arg  = 
            MC   = file_[0]
            prop = file_[1, 1:]
            analayzan(ax, filename, mean_E=boolians[i+3], mean_M=boolians[i+2],
                           spec_heat=boolians[i+1], susceptibility=boolians[i])
            ax.show()
            ax.clf()
    
    else:
        properties = np.array(["E", "M", "Cv", "X"])

        for element in properties:
            filename = np.loadtxt("{: s}.txt".format(element))

            analayzan(ax, filename)