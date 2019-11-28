all : solarsystem.out
	echo All done

test : visual_test_solver.out
	echo Test compilation done.

##############
# dependencies
##############

vector_3D.o : vector_3D.h vector_3D.cpp
	g++ -c vector_3D.cpp -std=c++17

solver.o : solver.h solver.cpp
	g++ -c solver.cpp -std=c++17

planet.o : planet.h planet.cpp
	g++ -c planet.cpp -std=c++17

#################
# programs to run
#################

solarsystem.o : planet.h solver.h vector_3D.h solver.o planet.o vector_3D.o
	g++ -c solarsystem.cpp -std=c++17

solarsystem.out : planet.h solver.h vector_3D.h solver.o planet.o vector_3D.o solarsystem.o
	g++ -o solarsystem.out solarsystem.o solver.o planet.o vector_3D.o -std=c++17

########
# tests
########

visual_test_solver.o : solver.h solver.o visual_test_solver.cpp
	g++ -c visual_test_solver.cpp -std=c++17

visual_test_solver.out : visual_test_solver.o solver.h solver.o
	g++ -o visual_test_solver.out visual_test_solver.o solver.o -std=c++17

##########
# clean up
##########

clean :
	rm *.o
	rm *.out
# 	rm solver.o planet.o solarsystem.o
# 	rm solarsystem.out

# clean_test :
# 	rm visual_test_solver.o
# 	rm visual_test_solver.out