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
    fetch_initial_parameters_from_file();

    // // U = {rx, ry, rz, vx, vy, vz}
    // // arma::vec U1 = {1, 0, 0, 0, 2*pi, 0};   // Earth
    // arma::vec U1 = fetch_initial_parameters_from_file();
    // arma::vec U2 = {0, 5.2, 0, 2.75, 0, 0}; // Jupiter
    // double mass1 = 5.972e24; // Earth mass.
    // double mass2 = 1.898e27; // Jupiter mass.
    // Solarsystem q;
    // q.add_planet(mass1, U1);
    // q.add_planet(mass2, U2);

    // // arma::vec u0 = q.get_U0();
    // // u0.print();
    // q.solve_system();
    // // arma::vec initial = q.get_U0();



    return 0;
}


