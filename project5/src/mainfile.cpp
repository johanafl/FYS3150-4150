#include "solarsystem.h"


arma::mat fetch_initial_parameters_from_file()
{   /*
    Fetch the initial parameters for every planet from a text file from
    NASA. Read values into an Armadillo matrix with planets as columns
    and parameters as rows.

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
    planet_data.open("init.txt", std::ios::in);

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


int main()
{   
    arma::mat all_planets_initial = fetch_initial_parameters_from_file();

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
    
    double sun_mass     = 1.9891e30;
    double mercury_mass = 3.285e23;
    double venus_mass   = 4.867e24;
    double earth_mass   = 5.972e24;
    double mars_mass    = 6.39e23;
    double jupiter_mass = 1.898e27;
    double saturn_mass  = 5.683e26;
    double uranus_mass  = 8.681e25;
    double neptune_mass = 1.024e26;
    double pluto_mass   = 1.309e22;
    
    Solarsystem q;
    // q.add_planet(sun_mass, sun_initial);
    // q.add_planet(mercury_mass, mercury_initial);
    // q.add_planet(venus_mass, venus_initial);
    q.add_planet(earth_mass, earth_initial);
    // q.add_planet(mars_mass, mars_initial);
    q.add_planet(jupiter_mass, jupiter_initial);
    // q.add_planet(saturn_mass, saturn_initial);
    // q.add_planet(uranus_mass, uranus_initial);
    // q.add_planet(neptune_mass, neptune_initial);
    // q.add_planet(pluto_mass, pluto_initial);

    q.solve_system();

    return 0;
}


