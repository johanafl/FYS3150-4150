#include "circular_matrix.h"

CircularMatrix::CircularMatrix(int n, double seed_input)
{   /*
    Parameters
    ----------
    n : int
        Dimension of matrix is n x n.
    
    seed_input : double
        Initializin the matrix radomly with given seed.
    */
    dim  = n;
    seed = seed_input;
    
    matrix = new double[dim*dim];
    initial_spin();
}

CircularMatrix::CircularMatrix(int n)
{   /*
    Parameters
    ----------
    n : int
        Dimension of matrix is n x n.
    
    seed_input : double
        Initializin the matrix radomly with Unix-time as seed.
    */
    time_t time_seed;
    time(&time_seed);
    
    seed = time_seed;
    dim  = n;

    matrix = new double[dim*dim];
    initial_spin();    
}

CircularMatrix::CircularMatrix(int n, double* init_set)
{   /*
    Parameters
    ----------
    n : int
        Dimension of matrix is n x n.
    
    init_set : double pointer
        Initializin the matrix with given input.
    
    Note
    ----
    Does not give error if input is larger or smaller than dimension n.
    */
    dim  = n;
    matrix = new double[dim*dim];

    for (int i = 0; i < dim*dim; i++)
    {   // populating the array with given array.
        matrix[i] = init_set[i];
    }
}

void CircularMatrix::initial_spin()
{   /*
    Initialize the matrix with spins drawn randomly using Mersenne
    Twister 19937 PRNG with seed from system time.
    */

    std::mt19937 engine(seed);
    std::uniform_int_distribution<int> uniform(0, 1);
    
    for (int i = 0; i < dim*dim; i++)
    {   // populating the array randomly with +- 1
        matrix[i] = 2*uniform(engine) - 1;
    }
}

void CircularMatrix::initial_spin(double seed_input)
{   /*
    Initialize the matrix with spins drawn randomly using Mersenne
    Twister 19937 PRNG with seed given as input.

    Parameters
    ----------
    seed_input : double
        Input seed value for the PRNG.
    */
    seed = seed_input;
    CircularMatrix::initial_spin();
}

void CircularMatrix::ordered_spin()
{   /*
    Initialize the spin matrix ordered with all spin up.
    */
    for (int i = 0; i < dim*dim; i++) {matrix[i] = 1;}
}

void CircularMatrix::set_new_dim_and_seed(int n, double new_seed)
{   /*
    Reset the spin matrix with size nxn and initial spins randomly
    drawn with input seed.

    Parameters
    ----------
    n : int
        Matrix dimension (nxn).

    new_seed : double
        New seed value for the PRNG.
    */
    delete[] matrix;
    dim  = n;
    seed = new_seed;
    
    matrix = new double[dim*dim];
    initial_spin();
}

void CircularMatrix::set_new_dim(int n)
{   /*
    Reset the spin matrix with size nxn.

    Parameters
    ----------
    n : int
        Matrix dimension (nxn).
    */
    CircularMatrix::set_new_dim_and_seed(n, seed);
}

void CircularMatrix::print()
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

double& CircularMatrix::operator() (int row, int col)
{   /*
    Index the matrix. Translate 2D indices to flat indices.

    Parameters
    ----------
    row : int
        Row index.

    col : int
        Column index.
    */
    
    return matrix[dim*((row + dim)%dim) + ((col + dim)%dim)];
}

double& CircularMatrix::operator() (int row, int col, bool safe)
{   /*
    Index the matrix. Translate 2D indices to flat indices. Boundary
    checks for matrix indexing. Allows index of padding, but nothing
    more.

    Parameters
    ----------
    row : int
        Row index.

    col : int
        Column index.

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

CircularMatrix::~CircularMatrix()
{
    delete[] matrix;
}
