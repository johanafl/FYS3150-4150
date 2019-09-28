all : test_jacobi.out quantum_dots.out quantum_dots_two_electrons.out
	
	echo All done

jacobi.o : jacobi.h jacobi.cpp

	g++ -c jacobi.cpp -std=c++17 -larmadillo

test_jacobi.o : jacobi.h test_jacobi.cpp

	g++ -c test_jacobi.cpp -std=c++17 -larmadillo

test_jacobi.out : jacobi.h jacobi.o test_jacobi.o

	g++ -o test_jacobi.out jacobi.o test_jacobi.o -std=c++17 -larmadillo

quantum_dots.o : jacobi.h quantum_dots.cpp

	g++ -c quantum_dots.cpp -std=c++17 -larmadillo

quantum_dots.out : jacobi.h quantum_dots.o jacobi.o

	g++ -o quantum_dots.out quantum_dots.o jacobi.o -std=c++17 -larmadillo 


quantum_dots_two_electrons.o : jacobi.h quantum_dots_two_electrons.cpp

	g++ -c quantum_dots_two_electrons.cpp -std=c++17 -larmadillo

quantum_dots_two_electrons.out : jacobi.h quantum_dots_two_electrons.o jacobi.o

	g++ -o quantum_dots_two_electrons.out quantum_dots_two_electrons.o jacobi.o -std=c++17 -larmadillo 

clean :

	rm test_jacobi.out quantum_dots.out quantum_dots_two_electrons.out
	rm test_jacobi.o quantum_dots.o jacobi.o quantum_dots_two_electrons.o