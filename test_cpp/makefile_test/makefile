# Define all the macros here

CC = g++
OBJ = jacobi.cpp
CFLAGS = -larmadillo

test.out : $(OBJ) jacobi.h test.cpp
	$(CC) $(OBJ) test.cpp $(CFLAGS) -o $@

run.out : $(OBJ) jacobi.h main.cpp
	$(CC) $(OBJ) main.cpp $(CFLAGS) -o $@

clean : 
	rm *.out