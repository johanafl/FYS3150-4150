all : energy_solver.out test_circular_matrix.out test_energy_solver.out
	echo All done

circular_matrix.o : circular_matrix.h circular_matrix.cpp
	g++ -c circular_matrix.cpp -std=c++17

test_circular_matrix.o : circular_matrix.h test_energy_solver.cpp
	g++ -c test_circular_matrix.cpp -std=c++17

test_circular_matrix.out : circular_matrix.h test_energy_solver.o circular_matrix.o
	g++ -o test_circular_matrix.out test_circular_matrix.o circular_matrix.o -std=c++17

test_energy_solver.o : circular_matrix.h test_energy_solver.cpp
	g++ -c test_energy_solver.cpp -std=c++17

test_energy_solver.out : circular_matrix.h test_energy_solver.o circular_matrix.o
	g++ -o test_energy_solver.out test_energy_solver.o circular_matrix.o -std=c++17

energy_solver.o : circular_matrix.h energy_solver.cpp
	g++ -c energy_solver.cpp -std=c++17

energy_solver.out : circular_matrix.h energy_solver.o circular_matrix.o
	g++ -o energy_solver.out energy_solver.o circular_matrix.o -std=c++17 -O3

clean :
	rm energy_solver.o circular_matrix.o test_circular_matrix.o test_energy_solver
	rm energy_solver.out test_circular_matrix.out test_energy_solver