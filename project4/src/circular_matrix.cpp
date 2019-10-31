#include <random>
#include <iostream>
#include <iomanip>
#include <stdlib.h>

class CircularMatrix
{
private:
    int dim;    // matrix is dim x dim
    int seed;   // seed for the random generator

public:
    double* matrix;
    
    CircularMatrix(int n, double seed_input)
    {   
        /*
        Parameters
        ----------
        n : int
            Dimension of matrix is n X n, so n i # rows (and/or # columns).
        
        seed_input : double
            Seed making the matrix random.
        */
        dim  = n;
        seed = seed_input;
        
        initial_spin();
    }

    CircularMatrix(int n)
    {   
        time_t time_seed;
        time(&time_seed);
        
        seed = time_seed;
        dim  = n;

        initial_spin();    
    }

    void initial_spin()
    {
        matrix = new double[dim*dim];
        std::mt19937 engine(seed);
        std::uniform_int_distribution<int> uniform(0, 1);
        
        for (int i = 0; i < dim*dim; i++)
        {   // populating the array randomly with +- 1
            matrix[i] = 2*uniform(engine) - 1;
        }
    }

    void print()
    {   /*
        Prints a nice visualization of the matrix.
        */ 
        std::cout << "[";
        
        for (int i = 0; i < dim*dim - 1; i++)
        {   // iterates over all elements in the matrix
            std::cout << std::setw(3) << matrix[i];
            if ( (i + 1)%dim == 0 )
            {   // prints a line shift for every new row
                std::cout << "\n" << " ";
            }
        }
        // special case for the last row
        std::cout << std::setw(3) << matrix[dim*dim-1] << "  ]" << std::endl;
    }

    double& operator() (int row, int col)
    {   /*
        Translates the 2D indices to flat indices.
        */
        return matrix[dim*(((row%dim) + dim)%dim) + (((col%dim) + dim)%dim)];
    }

    double& operator() (int row, int col, bool safe)
    {   /*
        Boundary checks for matrix indexing. Allows index of padding, but
        nothing more.
        */
        if ( (row < -1) or (row > dim) or (col < -1) or (col > dim) )
        {   // throws error if index is out of padding range  
            std::cout << "Indexing out of bounds. Terminating process." << std::endl;
            exit(EXIT_FAILURE);
        }
        else
        {   
            int row_idx = ((row%dim) + dim)%dim;
            int col_idx = ((col%dim) + dim)%dim;
            return matrix[dim*row_idx + col_idx];
        }
    }

    ~CircularMatrix()
    {
        delete[] matrix;
    }
};

// int main()
// {   
//     int seed = 1337;
//     int n = 3;
    
//     CircularMatrix q(n, seed);
//     q.print();
    
//     std::cout << q(0, 8, true) << std::endl;
//     return 0;
// }