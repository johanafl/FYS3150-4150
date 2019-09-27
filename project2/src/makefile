all : test.out quantum.out
	
	echo All done

jacobi.o : jacobi.h jacobi.cpp

	g++ -c jacobi.cpp -std=c++17 -larmadillo

test_jacobi.o : jacobi.h test_jacobi.cpp

	g++ -c test_jacobi.cpp -std=c++17 -larmadillo

test.out : jacobi.h jacobi.o test_jacobi.o

	g++ -o test.out jacobi.o test_jacobi.o -std=c++17 -larmadillo

quantum_dots.o : jacobi.h quantum_dots.cpp

	g++ -c quantum_dots.cpp -std=c++17 -larmadillo

quantum.out : jacobi.h quantum_dots.cpp jacobi.o

	g++ -o quantum.out quantum_dots.cpp jacobi.o -std=c++17 -larmadillo 

clean :

	rm test.out quantum.out 
	rm test_jacobi.o quantum.o jacobi.o