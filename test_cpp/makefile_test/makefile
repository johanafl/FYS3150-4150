# Really basic makefile syntax:
# The_name_of_the_file_that_is_produced(target/destination): the_file(s)_that_is_used_to_create_the_final_form(dependencies)
# 	The_shell_command_that_is_run_to_make_the_files(system command(s), needs to be one tab in)

my_output.exe: test_prog3.o
	g++ test_prog3.o -o my_output.exe

test_prog3.o: test_prog1.o test_prog2.o
	g++ -c test_prog3.cpp test_prog1.o test_prog2.o

test_prog1.o: test_prog1.cpp
	g++ -c test_prog1.cpp

test_prog2.o: test_prog1.cpp
	g++ -c test_prog1.cpp

# The target can also be an action that should be done, f.ex. 'clean':
# To run the clean, type in terminal$ make clean
clean:
	rm *.o