#include <random>
#include <iostream>
#include <iomanip>

class Matrix
{
private:
    int dim;    // matrix is dim x dim
    int seed;   // seed for the random generator


public:
    double *matrix;
    
    Matrix(int n, double seed_input)
    {   
        dim  = n;
        seed = seed_input;
        
        initial_spin();
    }

    Matrix(int n)
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

    void print_matrix()
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

    ~Matrix()
    {
        delete[] matrix;
    }
};



int main()
{   
    Matrix q(3);
    q.print_matrix();
    return 0;
}