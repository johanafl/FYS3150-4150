test.out: jacobi.h jacobi.cpp test_jacobi.cpp

	g++ -o test.out jacobi.cpp test_jacobi.cpp -std=c++17 -larmadillo

quantum.out: jacobi.h quantum_dots.cpp jacobi.cpp

	g++ -o quantum.out quantum_dots.cpp jacobi.cpp -std=c++17 -larmadillo -llapack -lblas

clean:

	rm test.out
	rm quantum.out