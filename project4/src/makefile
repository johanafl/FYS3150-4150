all : energy_solver.out
	echo All done

circular_matrix.o : circular_matrix.h circular_matrix.cpp
	g++ -c circular_matrix.cpp -std=c++17

energy_solver.o : circular_matrix.h energy_solver.cpp
	g++ -c energy_solver.cpp -std=c++17

energy_solver.out : circular_matrix.h energy_solver.o circular_matrix.o
	g++ -o energy_solver.out energy_solver.o circular_matrix.o -std=c++17

clean :
	rm energy_solver.o circular_matrix.o
	rm energy_solver.out