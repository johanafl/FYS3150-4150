# Project 4 - 

All program files are located in `src/`. There are six C++ files in total. The files circular_matrix.cpp and energy_solver.cpp will initialize the lattice and computing the different physical properties for the system. The file generate_parallel.cpp will parallelize the running of the program for a range of temperatures. The file generate_data.cpp will pluck your flower ;o) The two test files will test for different aspects of energy_solver.cpp and circular_matrix.cpp.

The Python file analyza_data.py will use the text files from the file that just plucked your flower to visualize the data.

Both the energy_solver.cpp and circular_matrix.cpp files be compiled the following way:
Linux:
`$  g++ program_name.cpp -o run.out -O3 -std=c++17`
macOS:
`$  clang++ program_name.cpp -o run.out -O3 -std=c++17`

The generate_data.cpp file with
generate_data:
`$  g++ -c generate_data.cpp -std=c++17; g++ -o generate_data.out generate_data.o energy_solver.o circular_matrix.o -std=c++17`

generate_data.o:
`$  g++ -c generate_data.cpp -std=c++17`

generate_data.out:
`$  g++ -o generate_data.out generate_data.o energy_solver.o circular_matrix.o -std=c++17`

The parallelized program uses (WE DIDN'T CHANGE ANYTHING HERE, DID WE?) MPIC2H, which can be downloaded at ´https://www.mpich.org/downloads/´. The program is compiled with
generate_parallel.o:
`$  mpic++ -c generate_parallel.cpp -std=c++17 -O3`

generate_parallel.out:
`$  mpic++ -o generate_parallel.out generate_parallel.o energy_solver.o circular_matrix.o -std=c++17 -O3`

The `doc/` directory contains the report for this project.
