#include "jacobi.h"
#include <fstream>
#include <iomanip>


arma::vec generate_diagonals(int n, double step, double rho_min)
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

    arma::vec diag(n);

    // rho = rho0 + ih
    // new: w*w*rho*rho + 1/rho

    // for (int i = 0; i < n; i++)
    // {
    //     diag(i) = w*w*(rho_min + i*step)*(rho_min + i*step) + 1/(rho_min + i*step);
    // }

    for (int i = 0; i < n; i++)
    {   // creating the diagonal elements
        diag(i) = 2/(step*step) + (rho_min + i*step)*(rho_min + i*step);
    }

    return diag;
}


class QuantumData
{

private:
    bool progress;      // boolean for toggling progress info on/off
    double eig;         // analytical eigenvalue
    double step;
    double off_diag;
    int n;                  // grid size
    double rho_min;
    double rho_max;         // approximating infinity
    double tol_off_diag;    // tolerance for Jacobi

    // loop-specific values
    double tmp;       // for reverting rho_max to original max value
    double d_rho;       // rho step size for incrementing in the loop
    double rho_end;        // end rho_max value for the loop
    int n_end;            // end grid value
    int dn;             // grid step size
    int num_eig;         // number of eigenvalues to write to file
    double num_rho; // number of rho_max values tested
    int num_n;
    std::ofstream data_file;
public:
    QuantumData()
    {
        progress = true;    // boolean for toggling progress info on/off
        n = 100;            // grid size
        rho_min = 0;
        rho_max = 5;         // approximating infinity
        tol_off_diag = std::pow(10, -5); // tolerance for Jacobi
        

        // loop-specific values
        tmp = rho_max;       // for reverting rho_max to original max value
        d_rho   = 0.1;       // rho step size for incrementing in the loop
        rho_end = 10;        // end rho_max value for the loop
        n_end = 410;            // end grid value
        dn = 10;             // grid step size
        num_eig = 8;         // number of eigenvalues to write to file
        num_rho = (rho_end - rho_max)/d_rho; // number of rho_max values tested
        num_n = (n_end - n)/dn;

        // generating data file
        std::ofstream data_file;
        data_file.open("eigenvalues2.txt", std::ios_base::app);
        data_file << num_eig << " " << num_rho << " " << num_n << "\n";
        data_file << std::setw(20) << "calculated";
        data_file << std::setw(20) << "exact";
        data_file << std::setw(20) << "error";
        data_file << std::setw(20) << "rho_max";
        data_file << std::setw(20) << "n\n";
    }

    void loop_rho()
    {   /*
        Loops over rho_max values.
        */
        
        while (rho_max < rho_end)
        {            
            eig = 3;    // initial analytical eigenvalue
            step = (rho_max - rho_min)/n;
            off_diag = -1/(step*step);

            // constructing vector for the diagonal elements            
            arma::vec diag = generate_diagonals(n, step, rho_min);

            // constructing tri-diagonal matrix
            arma::mat A = construct_diag_matrix(n, off_diag, diag);
            
            // using Jacobis method to extract the eigenvalues of the tri-diagonal matrix
            find_eig(n, A, tol_off_diag);
            arma::vec sorted_diag = arma::sort(A.diag(0));

            if (progress)
            {   // for printing progress data
                std::cout << "rho_max: " << rho_max << " rho_end: ";
                std::cout << rho_end << "  ";
                std::cout << "n: " << n << " n_end: ";
                std::cout << n_end << std::endl;
            }


            for (int i = 0; i < num_eig; i++)
            {   // writing data to file
                
                data_file << std::setw(20) << std::setprecision(10) << sorted_diag(i);
                data_file << std::setw(20) << std::setprecision(10) << eig;
                data_file << std::setw(20) << std::setprecision(10) << fabs(eig - sorted_diag(i));
                data_file << std::setw(20) << rho_max;
                data_file << std::setw(20) << n;
                data_file << "\n";
                eig = eig + 4;

            }
            rho_max = rho_max + d_rho;
        }
    }

    ~QuantumData()
    {
        data_file.close();
    }

};

void quantum_data(int task)
{   /*

    Loops over rho_max and n (grid) values, calculates eigenvalues by the Jacobi
    method, and writes data to a text file for visualization with Python.
    */



    bool progress = true;   // boolean for toggling progress info on/off
    double eig;             // analytical eigenvalue
    double step;
    double off_diag;
    int n = 100;                // grid size
    double rho_min = 0;
    double rho_max = 5;         // approximating infinity
    double tol_off_diag = std::pow(10, -5); // tolerance for Jacobi
    

    // loop-specific values
    double tmp = rho_max;       // for reverting rho_max to original max value
    double d_rho   = 0.1;       // rho step size for incrementing in the loop
    double rho_end = 10;        // end rho_max value for the loop
    int n_end = 410;            // end grid value
    int dn    = 10;             // grid step size
    int num_eig    = 8;         // number of eigenvalues to write to file
    double num_rho = (rho_end - rho_max)/d_rho; // number of rho_max values tested
    int num_n = (n_end - n)/dn;


    // generating data file
    std::ofstream data_file;
    data_file.open("eigenvalues2.txt", std::ios_base::app);
    data_file << num_eig << " " << num_rho << " " << num_n << "\n";
    data_file << std::setw(20) << "calculated";
    data_file << std::setw(20) << "exact";
    data_file << std::setw(20) << "error";
    data_file << std::setw(20) << "rho_max";
    data_file << std::setw(20) << "n\n";



    if (progress)
    {   // for printing progress data
        std::cout << "calculating" << std::endl;
        std::cout << "===========" << std::endl;
    }


    while (n < n_end)
    {   // looping over grid values
        
        while (rho_max < rho_end)
        {   // looping over rho_max values and writing data to file
            
            eig = 3;    // initial analytical eigenvalue
            step = (rho_max - rho_min)/n;
            off_diag = -1/(step*step);

    
            // constructing vector for the diagonal elements            
            arma::vec diag = generate_diagonals(n, step, rho_min);

            // constructing tri-diagonal matrix
            arma::mat A = construct_diag_matrix(n, off_diag, diag);
            
            // using Jacobis method to extract the eigenvalues of the tri-diagonal matrix
            find_eig(n, A, tol_off_diag);
            arma::vec sorted_diag = arma::sort(A.diag(0));

            if (progress)
            {   // for printing progress data
                std::cout << "rho_max: " << rho_max << " rho_end: ";
                std::cout << rho_end << "  ";
                std::cout << "n: " << n << " n_end: ";
                std::cout << n_end << std::endl;
            }


            for (int i = 0; i < num_eig; i++)
            {   // writing data to file
                
                data_file << std::setw(20) << std::setprecision(10) << sorted_diag(i);
                data_file << std::setw(20) << std::setprecision(10) << eig;
                data_file << std::setw(20) << std::setprecision(10) << fabs(eig - sorted_diag(i));
                data_file << std::setw(20) << rho_max;
                data_file << std::setw(20) << n;
                data_file << "\n";
                eig = eig + 4;

            }
            rho_max = rho_max + d_rho;
        }
        rho_max = 5;
        n += dn;
    }
    data_file.close();
}

int main()
{   

    QuantumData q;
    q.loop_rho();
    //quantum_data();
    return 0;
}