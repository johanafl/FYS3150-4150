#include "solver.h"
#include "planet.h"

double mass_earth   = 6e24;     // [kg]
double mass_jupiter = 1.9e27;   // [kg]
double mass_mars    = 6.6e23;   // [kg]
double mass_venus   = 4.9e24;   // [kg]
double mass_saturn  = 5.5e26;   // [kg]
double mass_mercury = 3.3e23;   // [kg]
double mass_uranus  = 8.8e25;   // [kg]
double mass_neptun  = 1.03e26;  // [kg]
double mass_pluto   = 1.31e22;  // [kg]

double dist_earth   = 1;       // [AU]
double dist_jupiter = 5.20;    // [AU]
double dist_mars    = 1.52;    // [AU]
double dist_venus   = 0.72;    // [AU]
double dist_saturn  = 9.54;    // [AU]
double dist_mercury = 0.39;    // [AU]
double dist_uranus  = 19.19;   // [AU]
double dist_neptun  = 30.06;   // [AU]
double dist_pluto   = 39.53;   // [AU]

void acceleration()
{
    // a = 0
    // for planet in planets:
    //      a += G*planet.mass/|planet.r - current_planet.r|^2
    // return a
}

int main()
{
    int num_planets = 2;
    Planet* planets = new Planet[num_planets];
    Planet earth(mass_earth, dist_earth, 0, 0, 0, 2*pi, 0);
    planets[0] = earth;
    
    // std::cout << earth.mass << std::endl;
    // std::cout << earth.init_x << std::endl;
    // std::cout << earth.init_y << std::endl;
    // std::cout << earth.init_z << std::endl;
    // std::cout << earth.init_vx << std::endl;
    // std::cout << earth.init_vy << std::endl;
    // std::cout << earth.init_vz << std::endl;
    Solver sovle_system(1e5, 1e-2, earth.init_x, earth.init_y, earth.init_vx, earth.init_vy);
    sovle_system.velocity_verlet();
    return 0;
}