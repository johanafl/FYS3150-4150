all : generate_data.out test_circular_matrix.out # energy_solver.out test_energy_solver.out
	echo All done

circular_matrix.o : circular_matrix.h circular_matrix.cpp
	g++ -c circular_matrix.cpp -std=c++17

test_circular_matrix.o : circular_matrix.h test_circular_matrix.cpp circular_matrix.o catch.hpp
	g++ -c test_circular_matrix.cpp -std=c++17

test_circular_matrix.out : circular_matrix.h test_circular_matrix.o circular_matrix.o catch.hpp
	g++ -o test_circular_matrix.out test_circular_matrix.o circular_matrix.o -std=c++17

# test_energy_solver.o : circular_matrix.h test_energy_solver.cpp
# 	g++ -c test_energy_solver.cpp -std=c++17 // REMEMBER CLEAN WHEN DONE!!!!

# test_energy_solver.out : circular_matrix.h test_energy_solver.o circular_matrix.o
# 	g++ -o test_energy_solver.out test_energy_solver.o circular_matrix.o -std=c++17

energy_solver.o : circular_matrix.h energy_solver.h energy_solver.cpp
	g++ -c energy_solver.cpp -std=c++17

# energy_solver.out : circular_matrix.h energy_solver.h energy_solver.o circular_matrix.o
# 	g++ -o energy_solver.out energy_solver.o circular_matrix.o -std=c++17 -O3

generate_data.o : circular_matrix.h energy_solver.h energy_solver.o circular_matrix.o
	g++ -c generate_data.cpp -std=c++17

generate_data.out : circular_matrix.h energy_solver.h energy_solver.o circular_matrix.o generate_data.o
	g++ -o generate_data.out generate_data.o energy_solver.o circular_matrix.o -std=c++17 -O3

clean :
	rm energy_solver.o circular_matrix.o generate_data.o test_circular_matrix.o
	rm generate_data.out test_circular_matrix.out
	# rm energy_solver.o circular_matrix.o test_circular_matrix.o test_energy_solver.o
	# rm energy_solver.out test_circular_matrix.out test_energy_solver