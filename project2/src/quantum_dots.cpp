#include "jacobi.h"
#include <fstream>
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
    */

    arma::vec diag(grid);

    for (int i = 0; i < grid; i++)
    {   // creating the diagonal elements
        // diag(i) = 2/(step*step) + (rho_min + i*step)*(rho_min + i*step);
        diag(i) = 2/(step*step) + freq*freq*(rho_min + i*step)*(rho_min + i*step) + coulomb/(rho_min + i*step);
    }

    return diag;
}


class QuantumData
{

private:
    bool progress;      // boolean for toggling progress info on/off
    bool looping_grid;
    bool looping_freq;
    double eig;         // analytical eigenvalue
    double step;
    double off_diag;
    int grid;              // grid size
    double rho_min;
    double rho_max;     // approximating infinity
    double tol_off_diag;// tolerance for Jacobi

    // loop-specific values
    double rho_tmp;         // for reverting rho_max to original max value
    double rho_end;     // end rho_max value for the loop
    int grid_tmp;          // for reverting n to original max value
    int grid_end;          // end grid value
    double d_rho;       // rho step size for incrementing in the loop
    int d_grid;             // grid step size
    int num_eig;        // number of eigenvalues to write to file per rho
    double num_rho;     // number of rho_max values tested
    int num_grid;          // number of grid values tested
    // std::ofstream data_file;

public:
    QuantumData()
    {
        progress = true;    // boolean for toggling progress info on/off
        looping_grid = false;
        looping_freq = false;
        grid = 100;            // grid size
        rho_min = std::pow(10, -7);
        rho_max = 5;        // approximating infinity
        tol_off_diag = std::pow(10, -5); // tolerance for Jacobi
        

        // loop-specific values
        rho_tmp = rho_max;   // for reverting rho_max to original max value
        grid_tmp = grid;     // for reverting rho_max to original max value
        d_rho   = 0.1;       // rho step size for incrementing in the loop
        d_grid = 10;         // grid step size
        rho_end = 5.6;       // end rho_max value for the loop
        grid_end = 130;      // end grid value
        num_eig = 8;         // number of eigenvalues to write to file
        num_rho = (rho_end - rho_max)/d_rho; // number of rho_max values tested
        num_grid = (grid_end - grid)/d_grid;

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
            find_eig(grid, A, tol_off_diag);
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

        looping_grid = true;

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
        looping_freq   = true;
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

    ~QuantumData()
    {
    }

};



int main()
{   

    QuantumData q;
    q.loop_frequency();
    //quantum_data();
    return 0;
}