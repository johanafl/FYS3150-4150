# Project 4 - 

All program files are located in `src/`. There are six C++ files in total. The files circular_matrix.cpp and energy_solver.cpp will initialize the lattice and compute the different physical properties for the system. The file generate_parallel.cpp will parallelize the running of the program for a range of temperatures. The file generate_data.cpp is the non-parallelized interface of the simulation (it is never used, just use generate_parallel.cpp instead). The two test files will test for different aspects of energy_solver.cpp and circular_matrix.cpp.

The Python file analyza_data.py uses the text files generated from the C++ programs to visually analyse the data.

Run the simulation by using the `makefile`. Compile by typing `make`, clean by typing `make clean`. Run the parallelized code by `mpiexec -n 12 ./generate_parallel.out` for running the program in parallel with 12 threads.

Remove the comments of the desired calculations in the main block of `analyse_data.py` and `generate_parallel.cpp`.

The `doc/` directory contains the report for this project as well as figures used in the report.
