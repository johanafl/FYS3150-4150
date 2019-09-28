#include "jacobi.h"
#include <iomanip>
#include <string>


arma::vec effective_potential(int grid, double step, double rho_min)
{
    /*
    Creates a vector containing correct diagonal elements for task 2d.

    Parameters
    ----------
    n : int
        Grid size.

    step : double
        Step size.

    rho_min : double
        Minimum rho value.
    */

    arma::vec diag(grid);

    for (int i = 0; i < grid; i++)
    {   // creating the diagonal elements
        diag(i) = 2/(step*step) + (rho_min + i*step)*(rho_min + i*step);
    }

    return diag;
}

// void exact()
// {
//     /* 
//     Exact radial wavefunction for two electrons. (NB! this is u(r), so not 
//     actually the wavefunction.)
    
//     Need r_vector (same as rho). Must be n long.
//     We calculate for l = 0, so set l = 0.
//     Need an n to loop over/set resolution.
//     */
//     double u[n];
//     // omega_r = 0.25 (corresponds to n = 2 in article by M. Taut.) -> eigenvalue (proportional to energy) lambda = 0.6250
//     for (i=0; i<n; i++)
//     {
//         r = r_vector[i]
//         u[i] = pow(r,l+1) * exp(-r*r/(8*(l + 1))) * (1 + r/(2*(l + 1)));
//     }

//     // omega_r = 0.05 (corresponds to n = 3 in article by M. Taut.) -> eigenvalue (proportional to energy) lambda = 0.1750
//     for (i=0; i<n; i++)
//     {
//         r = r_vector[i]
//         u[i] = pow(r,l+1) * exp(-r*r/(8*(4*l + 5))) * (1 + r/(2*(l + 1)) + r*r/(4*(l + 1)*(4*l + 5)));
//     }
// }

class QuantumData
{
    /* NEED TO COMMENT!*/
private:
    bool progress     = true;  // boolean for toggling progress info on/off
    bool looping_grid = false; // boolean for toggling progress of grid info on/off
    bool looping_freq = false; // boolean for toggling progress of frequency info on/off
    double eig;                // analytical eigenvalue
    double step;
    double off_diag;
    int grid       = 100;      // grid size
    double rho_min = std::pow(10, -7);
    double rho_max = 5;        // approximating infinity
    double tol_off_diag = std::pow(10, -10); // tolerance for Jacobi

    // loop-specific values
    double rho_tmp = rho_max; // for reverting rho_max to original max value
    double rho_end = 5.1;     // end rho_max value for the loop
    int grid_tmp   = grid;    // for reverting n to original max value
    int grid_end   = 120;     // end grid value
    double d_rho   = 0.1;     // rho step size for incrementing in the loop
    int d_grid     = 10;      // grid step size
    int num_eig    = 8;       // number of eigenvalues to write to file per rho
    double num_rho = (rho_end - rho_max)/d_rho; // number of rho_max values tested
    int num_grid   = (grid_end - grid)/d_grid;  // number of grid points tested

    std::ofstream data_file;

public:
    QuantumData()
    {
        // opening data file and writing title
        data_file.open("eigenvalues2.txt", std::ios_base::app);
        data_file << num_eig << " " << num_rho << " " << num_grid << "\n";
        data_file << std::setw(20) << "calculated";
        data_file << std::setw(20) << "exact";
        data_file << std::setw(20) << "error";
        data_file << std::setw(20) << "rho_max";
        data_file << std::setw(20) << "n\n";
    }

    void set_grid_values(int grid_input, int grid_end_input, double d_grid_input)
    {   /* 
        For setting grid values. If not called, values are set to default.

        Parameters
        ----------
        grid_input : int
            Grid point start value.

        grid_end_input : int
            Grid point end value.

        d_grid_input : double
            Grid point step size
        */
        
        grid = grid_input;
        grid_end = grid_end_input;
        d_grid = d_grid_input;
        grid_tmp = grid;
        num_grid = (grid_end - grid)/d_grid;
    }

    void set_rho_values(double rho_max_input, double rho_end_input, double d_rho_input)
    {   /*

        For setting rho values. If not called, values are set to default.
        
        Parameters
        ----------
        rho_max_input : double
            rho max start value.

        rho_end_input : double
            rho max end value.

        d_rho_input : double
            rho max step size.
        */
        rho_max = rho_max_input;
        rho_tmp = rho_max;
        rho_end = rho_end_input;
        d_rho   = d_rho_input;
        num_rho = (rho_end - rho_max)/d_rho;
    }

    void set_progress_values(bool progress_input)
    {
        /* REMEMBER TO COMMENT*/
        progress = progress_input;
    }

    void loop_rho()
    {   /*
        Loops over rho_max values.
        */
       
        if ((progress) && (!looping_grid))
        {   // progress information
           std::cout << "looping over rho max" << std::endl;
        }
        
        while (rho_max < rho_end)
        {            
            eig = 3;    // initial analytical eigenvalue
            step = (rho_max - rho_min)/grid;
            off_diag = -1/(step*step);

            // constructing vector for the diagonal elements            
            arma::vec diag = effective_potential(grid, step, rho_min);

            // constructing tri-diagonal matrix
            arma::mat A = construct_diag_matrix(grid, off_diag, diag);
            
            // using Jacobis method to extract the eigenvalues of the tri-diagonal matrix
            arma::mat R = find_eig(grid, A, tol_off_diag);
            arma::vec sorted_diag = arma::sort(A.diag(0));

            if (progress)
            {   // for printing progress data
                std::cout << "rho_max: " << rho_max << " rho_end: ";
                std::cout << rho_end << "  ";
                std::cout << "n: " << grid << " n_end: ";
                std::cout << grid_end << std::endl;
            }

            for (int i = 0; i < num_eig; i++)
            {   // writing data to file
                data_file << std::setw(20) << std::setprecision(10) << sorted_diag(i);
                data_file << std::setw(20) << std::setprecision(10) << eig;
                data_file << std::setw(20) << std::setprecision(10) << fabs(eig - sorted_diag(i));
                data_file << std::setw(20) << rho_max;
                data_file << std::setw(20) << grid;
                data_file << "\n";
                eig = eig + 4;
            }
            rho_max = rho_max + d_rho;
        }
    }

    void loop_grid()
    {   /*
        Loops over grid (n) values.
        */

        if ((progress) && (!looping_freq))
        {   // progress information
            std::cout << "looping over grid and rho max" << std::endl;
        }

        while (grid < grid_end)
        {
            loop_rho();
            
            rho_max = rho_tmp;
            grid   += d_grid;
        }
    }

    ~QuantumData()
    {
        data_file.close();
    }

};

int main()
{   
    QuantumData q;
    q.set_grid_values(100, 101, 10);
    q.set_rho_values(1, 7, 0.1);
    q.loop_grid();

    return 0;
}