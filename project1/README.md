# Project 1 - Introduction to C++ and linear algebra problems in Computational Physics

We begin the course by studying how we can solve the Poisson differential equation using matricies and vector algebra. Spesificly, we will use the Thomas algorithm, a spesialized verison of the Thomas algorithm and LU decomposition to solve the problem. Since all of the above algorithms can be quite computationaly heavy, we will use the programming language C++ to to do the actual calculation and Python to visualize it. 

We will also write a report on this exersice. This report will be written as if it was not an exersice and instead as an actual scientific report.

The code consists of two files; algorithm_calculations.cpp and algorithm_analysis.py. To run the code, first compile the program with `g++ algorithm_calculations.cpp -larmadillo -std=c++17 -o run.out` and run it with `./run.out`. All data is calculated and text files needed for visualization are created. Run then the python file `python3 algorithm_analysis.py` to show all plots.

macOS users may run into a unresolved armadillo error when running `calculate_error();` in the C++ file. The error may look like this:

```bash
warning: solve(): system seems singular; attempting approx solution
** On entry to DLASCL, parameter number  4 had an illegal value
** On entry to DLASCL, parameter number  4 had an illegal value
'./run.out' terminated by signal SIGSEGV (Address boundary error)

warning: solve(): system seems singular (rcond: 6.94233e-79); attempting approx solution
```

We don't have a solution for this at the moment, other than trying to run `calculate_error();` 3-4 times until it shows no errors (it is somewhat random if it works).