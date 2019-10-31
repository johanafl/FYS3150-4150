all : energy_solver.out
	echo All done

energy_solver.o : energy_solver.cpp
	g++ -c energy_solver.cpp -std=c++17

circular_matrix.o : circular_matrix.cpp
	g++ -c circular_matrix.cpp -std=c++17

energy_solver.out : energy_solver.o circular_matrix.o
	g++ -o energy_solver.out energy_solver.o circular_matrix.o -std=c++17

clean :
	rm energy_solver.o circular_matrix.o
	rm energy_solver.out