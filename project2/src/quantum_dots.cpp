#include "jacobi.h"
#include <iomanip>
#include <string>


arma::vec generate_diagonals(int grid, double step, double rho_min, double freq, int coulomb)
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

    freq : double
        A freqency in front of the harmonic oscillator potential. 
        For the case of two electrons this normalizes the constants in the 
        coulomb potential and makes it unitless.
    
    coulomb : int, either 1 or 0
        A variable to indicate if we want to use an effictive potential with the
        coulomb potential. 
        Set to 1 to include.
        Set to 0 to exclude.
    */

    arma::vec diag(grid);

    for (int i = 0; i < grid; i++)
    {   // creating the diagonal elements
        diag(i) =   2/(step*step) 
                  + freq*freq*(rho_min + i*step)*(rho_min + i*step) // Harmonic oscillator potential (with frequency)
                  + coulomb/(rho_min + i*step);                     // Coulomb potential
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
//     // omega_r = 4 (corresponds to n = 2 in article by M. Taut.) -> eigenvalue (proportional to energy) lambda = 0.6250
//     for (i=0; i<n; i++)
//     {
//         r = r_vector[i]
//         u[i] = pow(r,l+1) * exp(-r*r/(8*(l + 1))) * (1 + r/(2*(l + 1)));
//     }

//     // omega_r = 20 (corresponds to n = 3 in article by M. Taut.) -> eigenvalue (proportional to energy) lambda = 0.1750
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
    double tol_off_diag = std::pow(10, -5); // tolerance for Jacobi

    // loop-specific values
    double rho_tmp = rho_max; // for reverting rho_max to original max value
    double rho_end = 5.6;     // end rho_max value for the loop
    int grid_tmp   = grid;    // for reverting n to original max value
    int grid_end   = 130;     // end grid value
    double d_rho   = 0.1;     // rho step size for incrementing in the loop
    int d_grid     = 10;      // grid step size
    int num_eig    = 8;       // number of eigenvalues to write to file per rho
    double num_rho = (rho_end - rho_max)/d_rho; // number of rho_max values tested
    int num_grid   = (grid_end - grid)/d_grid;  // number of grid points tested

public:
    QuantumData()
    {
        /* Empty right now. Might include accesibility for useres later.*/
    }

    void set_grid_values(int n, int end, double stepsize)
    {
        /* REMEMBER TO COMMENT*/
        grid = n;
        grid_end = end;
        d_grid = stepsize;
        grid_tmp = grid;
        num_grid = (grid_end - grid)/d_grid;
    }

    void set_rho_values(double my_man_rho, double the_rho_has_landed, double see_you_later_rho_i_gator, double sthep_that_steip)
    {
        /* REMEMBER TO COMMENT*/
        rho_min = my_man_rho;
        rho_max = the_rho_has_landed;
        rho_tmp = rho_max;
        rho_end = see_you_later_rho_i_gator;
        d_rho   = sthep_that_steip;
        num_rho = (rho_end - rho_max)/d_rho;
    }

    void set_progress_values(bool do_you_really_want_to_see_this, bool of_course_i_do, bool or_mayby_not)
    {
        /* REMEMBER TO COMMENT*/
        progress     = do_you_really_want_to_see_this;
        looping_grid = of_course_i_do;
        looping_freq = or_mayby_not;
    }

    void write_title_to_file(std::ofstream &data_file)
    {
        // opening data file and writing title
        // data_file.open("eigenvalues2.txt", std::ios_base::app);
        data_file << num_eig << " " << num_rho << " " << num_grid << "\n";
        data_file << std::setw(20) << "calculated";
        data_file << std::setw(20) << "exact";
        data_file << std::setw(20) << "error";
        data_file << std::setw(20) << "rho_max";
        data_file << std::setw(20) << "n\n";
    }

    void loop_rho(std::ofstream &data_file, double freq, int coulomb)
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
            arma::vec diag = generate_diagonals(grid, step, rho_min, freq, coulomb);

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

    void loop_grid(std::ofstream &data_file, double freq, int coulomb)
    {   /*
        Loops over grid (n) values.
        */

        if ((progress) && (!looping_freq))
        {   // progress information
            std::cout << "looping over grid and rho max" << std::endl;
        }

        while (grid < grid_end)
        {
            loop_rho(data_file, freq, coulomb);
            
            rho_max = rho_tmp;
            grid   += d_grid;
        }
    }

    void loop_grid()
    {   /*
        Loops over grid (n) values.
        */

        std::ofstream data_file;
        data_file.open("eigenvalues2.txt", std::ios_base::app);
        loop_grid(data_file, 1, 0);
        data_file.close();
    }

    void loop_frequency()
    {   /*
        Loops over frequencies (w).
        */
    
        std::string filename;
        std::ofstream data_file;
        double freq[4] = {0.01, 0.5, 1, 5};

        if (progress)
        {   // progress information
            std::cout << "looping over frequency, grid and rho max" << std::endl;
        }
 
        for (int i = 0; i < 4; i++)
        {   
            filename = "eigenvalues_w_" + std::to_string(freq[i]) + ".txt";
            data_file.open(filename, std::ios_base::app);
            write_title_to_file(data_file);
            
            loop_grid(data_file, freq[i], 1);
            
            grid = grid_tmp;
            data_file.close();
        }

    }

};

int main()
{   
    QuantumData q;
    q.loop_frequency();

    return 0;
}