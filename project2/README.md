# Project 2 - Solving eigenvalue problems to calculate the radial wavefunction of one and two electrons in a sperically symmetric harmonic oscillator potential.

All program files are located in `src/`. There are two C++ files which generate all the data needed for the analysis. The file `quantum_dots.cpp` generate all the data for the single electron system and writes the data to the textfile `eigenvalues.txt`. The second file `quantum_dots_two_electrons.cpp` generates all the data for the two electron system and writes the data to a set of files, `eigenvector_omega_<frequency>.txt` and `eigenvalue_omega_<frequency>.txt`. Run the C++ files by using the `make` file present in the same directory. The make file creates three executables, `quantum_dots.out`, `quantum_dots_two_electrons.out`, `test_jacobi.out`. The test checks for conservation of orthogonality, tests that the correct max value is found, and tests that the inner product is conserved. To generate all data, run:

```
make

./quantum_dots.out
./quantum_dots_two_electrons.out
```

The single Python file, `quantum_dots.py` does the analysis for both systems. The class `VisualizeData` does all the visualization for the single electron system, and all the function calls to generate plots are located in the main block at the bottom of the file. The functions `visualize_eigendata_two_electrons_numerical_and_analytical()` and `visualize_eigendata_two_electrons_numerical()` does the number crunching for the two electron system. Run the Python file by

```
python3 quantum_dots.py
```


fin