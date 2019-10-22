# Project 3 - Integrating quantumechanical integral with Gaussian quadrature and Monte Carlo integration

All program files are located in `src/`. There are six C++ files in total. The file gauss_laguerre.cpp is supplied by Morten to help us compute the wheights of for the Gaussian quadrature. The files gauss_legendre.cpp and gauss_laguerre.cpp evaluates the integral using Gaussian quadrature with their respective orthogonal polynomial basis. The files mc_*.cpp files evaluates the integral using Monte Carlo integration. 

All cpp-files, exept gauss_laguerre.cpp, create text files, which are placed int a separate directory, data_files. The Python file analyze_data.py uses these text files to visualize the data.

All the .cpp files save for the parallelized Monte Carlo can be compiled the following way:
Linux:
`$  g++ program_name.cpp -o run.out -O3 -std=c++17`
macOS:
`$  clang++ program_name.cpp -o run.out -O3 -std=c++17`

The parallelized program uses MPIC2H, which can be downloaded at ´https://www.mpich.org/downloads/´. The program is compiled with
`mpic++ program_name.cpp -o run.out -O3 -std=c++17`
and to run it with 10 threads,
`mpiexec -n 10 run.out`