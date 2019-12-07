all : main.out

tests : test_solver.out test_solar_system.out

test_solar_system.out : solver.h solar_system.h test_solar_system.cpp
	g++ -o test_solar_system.out test_solar_system.cpp -std=c++17

tests_solver.out : solver.h test_solver.cpp
	g++ -o test_solver.out test_solver.cpp -std=c++17

main.out : solver.h solar_system.h main.cpp
	g++ -o main.out main.cpp -std=c++17
	
clean :
	rm *.out
