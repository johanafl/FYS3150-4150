# Project 5 - Solar system

## Program structure
### solver.h
The solver is the bedrock of this project and it is modelled after ```ODESolver.py``` from the course IN1900 (previously INF1100). The layout is as follows: A superclass ```Solver``` contains general methods needed for solving ordinary differential equations. Implementations of specific algorithms are made as subclasses which inherit from ```Solver```. The subclasses contain only an advance method which advances the integration by one step, using a specific ODE solving algorithm. Say that we want to implement a Runge Kutta method of fourth order. We then simply need to make a class ```RungeKutta4``` which contain a single advance method with the Runge Kutta 4 algorithm. The Runge Kutta class inherits everything else it needs from the ```Solver``` superclass.

Advantages: Little coding is needed to implement a specific algorithm, since the framework is already present in the ```Solver``` class. It is elegant and easy to read.

Disadvantges: It is difficult and time consuming to write something as general as this. A system specifically designed for solving the solar system project would definitely be faster to implement, since we can taylor everything specifically to the project.

The contents of the ```Solver``` class are as follows: The constructor takes as input the number of integration steps and the number of stellar objects. (The final goal for us was to create something 100% general, so taking the number of stellar objects as input breaks this generality, but it is present simply because we didn't have enough time to generalise everything.) A method ```set_initial_conditions``` takes an ```arma::vec``` as input, which contains all the initial positions and initial velocities. The layout is vec = {rx1, ry1, rz1, vx1, vy1, vz1, ..., rxN, ryN, rzN, vxN, vyN, vzN}, where N is the number of stellar objects. ```Solver``` then contains a method ```solve``` which loops over the number of integration steps and calls a method ```advance``` for every iteration. ```Solver```'s ```advance``` method will be overwritten by the subclass' ```advance``` method, and it must therefore be ```virtual```. A method ```write_to_file``` writes all the data to a text file, where a column represents each planets x, y, z position and velocity with the same layout as the initial conditions vector. Another method ```write_selection_to_file``` writes only the data from the zeroth step to ```selection_1```th step, and the data from the ```selection_2```th step to the end. In other words, it writes only a selection at the beginning of the data set and a selection at the end of the data set. This is for the simulations which have a long simulation time and a low step length, where we are only interested in a portion of the data.

A separate implementation of the Sun-Mercury system with the relatiistic correction term had to be implemented as we havent had time to generalise this yet. It consists equally of a ```solve_mercury``` method.

The contents of the subclasses ```ForwardEuler``` and ```VelocityVerlet``` need only to contain ```advance``` methods for the specific algorithms, since ```Solver``` deals with the rest. The ```advance``` methods solve the acceleration, velocity and position for the next step. The actual function which is integrated (the RHS of the ODE, called ```acceleration```) is passed to the ```advance``` methods as an argument contained in a class object. The original plan was to pass the function as a pointer, but since the ```acceleration```s are member functions, we can't pass it as pointers, so we resolved the problem by calling by reference for the class object of which ```acceleration``` is contained.

### solar_system.h
The ```solar_system.h``` file contains all the specifics for keeping track on the solar system and for solving the solar system. The class ```SolarSystem``` has public methods ```add_celestial_body``` and ```solve_system```. ```add_celestial_body``` is for adding celestial bodies (planets, moons, stars) to the solar system. It takes as input the mass of the object as a ```double``` and the initial x, y, z position and velocity as an ```arma::vec```. The mass, position and velocity initial inputs are added to three internal arrays for storing all initial conditions for all added celestial bodies. The class ```SolarSystem``` also contains a simple ```resize``` method for resizing the initial condition arrays should they fill up. When all celestial bodies are added, the user enables the ```solve_system``` public method which in turn does a bunch of stuff: It initialises a class object of the desired ODE solving algorithm (Velocity Verlet or Forward Euler). It then fetches the initial conditions for all celestial bodies from the ```SolarSystem``` class using the private ```get_U0``` method, and inserts this information into the ```set_initial_conditions``` method of the ```Solver``` class. Then, the ```solve``` method of the ```Solver``` class is called which solves the ODE for the given amount of steps. Then, the ```write_to_file``` method may be called if we wish to save the data.

A specific example is:

```c++
SolarSystem q;
q.add_celestial_body(sun_mass, sun_initial);
q.add_celestial_body(mercury_mass, mercury_initial);
```
for initialising the ```SolarSystem``` class and adding celestial bodies, and then 

```c++
q.solve_system(num_steps, dt, func_id, method, outfilepath); // Solves and writes data to file.

q.solve_system(num_steps, dt, func_id, method); // Only solving, no writing.
```