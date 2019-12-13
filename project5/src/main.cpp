#include "solar_system.h"

// All masses in kg.
const double sun_mass     = 1.9891e30;
const double mercury_mass = 3.285e23;
const double venus_mass   = 4.867e24;
const double earth_mass   = 5.972e24;
const double mars_mass    = 6.39e23;
const double jupiter_mass = 1.898e27;
const double saturn_mass  = 5.683e26;
const double uranus_mass  = 8.681e25;
const double neptune_mass = 1.024e26;
const double pluto_mass   = 1.309e22;

arma::mat fetch_initial_parameters_from_file(std::string filepath)
{   /*
    Fetch the initial parameters for every planet from a text file from
    NASA. Read values into an Armadillo matrix with planets as columns
    and parameters as rows.

    Parameters
    ----------
    filepath : std::string
        The path and name of the text file.

    Returns
    -------
    stellar_init : arma::mat
        Matrix of initial parameters for all planets (and Pluto).
    */

    int num_input_params = 6;   // Number of input paramerers per stellar object.

    std::string line;
    std::string word;
    arma::mat stellar_init(num_input_params, 10);
    std::ifstream planet_data;
    planet_data.open(filepath, std::ios::in);

    // First indices of every the number and the length of every number.
    int number_length = 20;
    arma::Col<int> indices = {12, 37, 62, 87, 112, 137};
    
    for (int line_number = 0; line_number < 10; line_number++)
    {   // Iterate over each line in the file.
        
        std::getline(planet_data, line);
        
        for (int param = 0; param < num_input_params; param++)
        {   // Loop over rx, ry, rz, vx, vy, vz (input params).
            
            for (int i = indices(param); i <= indices(param) + number_length; i++)
            {   // Loop over all characters in each number.
                
                word += line[i];
            }
            
            // Convert the word to double and insert into matrix.
            stellar_init(param, line_number) = std::stod(word);
            word.clear();
        }
    }
    planet_data.close();

    // Convert the velocities from AU/d to AU/yr.
    stellar_init.submat(3, 0, 5, 9) *= 365.25;

    return stellar_init;
}


void task_5c()
{   /*
    Eath with Sun fixed at the origin.
    */
    arma::vec earth_initial = {1, 0, 0, 0, 2*pi, 0};

    double dt[4] = {1e-3, 1e-2};
    int func_id = 2;

    for (int i = 0; i < 2; i++)
    {   
        int num_steps = 100/dt[i];
        std::cout << "Generating task 5c data. dt = " << std::to_string(dt[i]) << std::endl;
        std::string filepath_fe = "data_files/task_5c_fe_dt=" + std::to_string(dt[i]) + ".txt";
        std::string filepath_vv = "data_files/task_5c_vv_dt=" + std::to_string(dt[i]) + ".txt";
        std::string method_fe = "Forward Euler";
        std::string method_vv = "Velocity Verlet";
        
        SolarSystem q;
        q.add_celestial_body(earth_mass, earth_initial);
        
        q.solve_system(num_steps, dt[i], func_id, method_fe, filepath_fe);
        q.solve_system(num_steps, dt[i], func_id, method_vv, filepath_vv);

        std::cout << std::endl;
    }
}


void task_5c_algorithm_timing()
{   
    int func_id = 2;    // Choose which of the accelerations to use.
    double dt = 1e-3;
    double simulation_time_in_years = 100000;
    int num_steps = simulation_time_in_years/dt;
    arma::vec earth_initial = {1, 0, 0, 0, 2*pi, 0};

    SolarSystem q;
    q.add_celestial_body(earth_mass, earth_initial);
    
    std::string method = "Forward Euler";
    std::cout << "Method: " << method << ". dt: " << dt << " yr, N: "
    << num_steps << ", T: " << simulation_time_in_years << " yr.\n" << std::endl;
    q.solve_system(num_steps, dt, func_id, method);

    method = "Velocity Verlet";
    std::cout << "Method: " << method << ". dt: " << dt << " yr, N: "
    << num_steps << ", T: " << simulation_time_in_years << " yr.\n" << std::endl;
    q.solve_system(num_steps, dt, func_id, method);
}


void task_5d()
{   /*
    Earth with Sun fixed at the origin. Earth with various initial
    velocities.
    */

    int func_id = 2;
    double dt = 1e-4; // Irrelevant since we only consider the initial energy.
    int num_steps = 1;
    arma::vec earth_initial = {1, 0, 0, 0, 0, 0};
    std::string method = "Velocity Verlet";
    std::string filepath;
    std::string tmp;
    
    int file_counter = 0;
    
    for (double i = 1; i <= 1.6; i = i + 0.00625)
    {   
        earth_initial[4] = i*2*pi;
        SolarSystem q;
        filepath = "data_files/task_5d_" + std::to_string(earth_initial[4]) + ".txt";
        
        q.add_celestial_body(earth_mass, earth_initial);
        q.solve_system(num_steps, dt, func_id, method, filepath);

        file_counter++;
    }
}


void task_5d_beta()
{   /*
    Earth with Sun fixed at the origin. Varying beta.
    */

    int func_id = 3;
    double dt = 1e-5;
    double simulation_time_in_years = 40;
    int num_steps = simulation_time_in_years/dt;
    std::cout << "task_5d_beta" << std::endl;
    std::cout << "dt: " << dt << " yr, simulation time: "
    << simulation_time_in_years << " yr, steps: " << num_steps << std::endl;
    
    arma::vec earth_initial = {1, 0, 0, 0, 6.7, 0};
    std::string method = "Velocity Verlet";
    std::string filepath;

    SolarSystem q;
    q.add_celestial_body(earth_mass, earth_initial);

    double beta[4] = {2, 2.33333333, 2.66666667, 3};
    
    for (int i = 0; i < 4; i++)
    {   
        filepath = "data_files/varying_beta=" + std::to_string(beta[i]) + ".txt";
        
        q.set_beta(beta[i]);
        q.solve_system(num_steps, dt, func_id, method, filepath);

    }
}


void task_5e()
{   /*
    Lets see how much Jupiter affects Earths orbit. Keeping the Sun
    fixed as the CM, using acceleration_2.
    */

    std::cout << "task_5e\n" << std::endl;
    
    std::string infilepath = "data_files/initial_parameters_solar_system_sun_at_rest.txt";
    arma::mat all_planets_initial = fetch_initial_parameters_from_file(infilepath);

    arma::vec earth_initial   = all_planets_initial.col(3);
    arma::vec jupiter_initial = all_planets_initial.col(5);
    
    int func_id = 2;    // Acceleration with the Sun fixed at the center.
    double dt = 1e-5;
    double simulation_time_in_years = 30;
    int num_steps = simulation_time_in_years/dt;
    // int num_steps = 10;
    
    std::string method = "Velocity Verlet";
    std::string filepath;

    double jupiter_factor[4] = {1, 10, 100, 1000};

    for (int i = 0; i < 4; i++)
    {
        SolarSystem q;
        q.add_celestial_body(earth_mass, earth_initial);
        q.add_celestial_body(jupiter_mass*jupiter_factor[i], jupiter_initial);
        
        filepath = "data_files/jupiter_mass=" + std::to_string(jupiter_mass*jupiter_factor[i]) + ".txt";
        q.solve_system(num_steps, dt, func_id, method, filepath);
    }

}


void task_5f()
{
    
    std::string infilepath = "data_files/initial_parameters_solar_system_sun_at_rest.txt";
    arma::mat all_planets_initial = fetch_initial_parameters_from_file(infilepath);

    arma::vec sun_initial     = all_planets_initial.col(0);
    arma::vec mercury_initial = all_planets_initial.col(1);
    arma::vec venus_initial   = all_planets_initial.col(2);
    arma::vec earth_initial   = all_planets_initial.col(3);
    arma::vec mars_initial    = all_planets_initial.col(4);
    arma::vec jupiter_initial = all_planets_initial.col(5);
    arma::vec saturn_initial  = all_planets_initial.col(6);
    arma::vec uranus_initial  = all_planets_initial.col(7);
    arma::vec neptune_initial = all_planets_initial.col(8);
    arma::vec pluto_initial   = all_planets_initial.col(9);
    
    const int func_id = 4;    // Acceleration with the Sun fixed at the center.
    const double dt = 1e-5;
    const double simulation_time_in_years = 100;
    const int num_steps = simulation_time_in_years/dt;

    std::cout << "task_5f" << std::endl;
    std::cout << "dt: " << dt << " yr, simulation time: "
    << simulation_time_in_years << " yr, steps: " << num_steps << std::endl;
    
    std::string method = "Velocity Verlet";
    std::string filepath = "data_files/sun_earth_jupiter.txt";

    SolarSystem q;
    q.add_celestial_body(sun_mass, sun_initial);
    q.add_celestial_body(earth_mass, earth_initial);
    q.add_celestial_body(jupiter_mass, jupiter_initial);
    q.solve_system(num_steps, dt, func_id, method, filepath);
}


void task_5g()
{   
    double dt = 1e-5;
    int num_steps = 2/dt;
    // arma::vec mercury_initial = {-3.397162482844107e-1, 1.179729441938765e-1, 4.080404839815078e-2, -1.504265080337234e-2*365.242199, -2.537738274969423e-2*365.242199, -6.937544462359349e-4*365.242199};
    arma::vec mercury_initial = {0.3075, 0, 0, 0, 12.44, 0};
    // arma::vec mercury_initial = {0.3075, 0, 0, 0, 13.9, 0};

    std::string filepath = "data_files/task_5g.txt";
    
    SolarSystem q;
    q.add_celestial_body(mercury_mass, mercury_initial);
    q.solve_system_mercury(num_steps, dt, filepath);
}


void task_5f_all_planets()
{   
    const double dt = 1e-3;
    const double simulation_time_in_years = 50;
    const int num_steps = simulation_time_in_years/dt;
    const int func_id = 4;        // Acceleration with cm at rest.
    std::string infilepath = "data_files/initial_parameters_solar_system_cm_at_rest.txt";
    arma::mat all_planets_initial = fetch_initial_parameters_from_file(infilepath);

    std::cout << "task_5f_all_planets" << std::endl;
    std::cout << "dt: " << dt << " yr, simulation time: "
    << simulation_time_in_years << " yr, steps: " << num_steps << std::endl;

    // All positions in AU, velocities in AU/yr.
    const arma::vec sun_initial     = all_planets_initial.col(0);
    const arma::vec mercury_initial = all_planets_initial.col(1);
    const arma::vec venus_initial   = all_planets_initial.col(2);
    const arma::vec earth_initial   = all_planets_initial.col(3);
    const arma::vec mars_initial    = all_planets_initial.col(4);
    const arma::vec jupiter_initial = all_planets_initial.col(5);
    const arma::vec saturn_initial  = all_planets_initial.col(6);
    const arma::vec uranus_initial  = all_planets_initial.col(7);
    const arma::vec neptune_initial = all_planets_initial.col(8);
    const arma::vec pluto_initial   = all_planets_initial.col(9);

    
    SolarSystem q;
    q.add_celestial_body(sun_mass, sun_initial);
    q.add_celestial_body(mercury_mass, mercury_initial);
    q.add_celestial_body(venus_mass, venus_initial);
    q.add_celestial_body(earth_mass, earth_initial);
    q.add_celestial_body(mars_mass, mars_initial);
    q.add_celestial_body(jupiter_mass, jupiter_initial);
    q.add_celestial_body(saturn_mass, saturn_initial);
    q.add_celestial_body(uranus_mass, uranus_initial);
    q.add_celestial_body(neptune_mass, neptune_initial);
    q.add_celestial_body(pluto_mass, pluto_initial);


    std::string outfilepath = "data_files/all_planets.txt";
    std::string method = "Velocity Verlet";
    q.solve_system(num_steps, dt, func_id, method, outfilepath);
}


int main()
{   
    // task_5c();
    // task_5c_algorithm_timing();
    // task_5d();
    // task_5d_beta();
    // task_5e();
    // task_5f();
    // task_5f_all_planets();
    task_5g();


    return 0;
}


