all : mainfile.out
	echo All done

test : visual_test_solver.out
	echo Test compilation done.

##############
# dependencies
##############

solver.o : solver.h solver.cpp
	g++ -c solver.cpp -std=c++17

#################
# programs to run
#################

solarsystem.o : solver.h solver.o solarsystem.h
	g++ -c solarsystem.cpp -std=c++17

mainfile.out : mainfile.cpp solver.h solver.o solarsystem.h solarsystem.o
	g++ -o mainfile.out mainfile.cpp solarsystem.o solver.o -std=c++17

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
