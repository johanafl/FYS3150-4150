#include "solarsystem.h"

int main()
{   
    arma::vec U0 = {1, 0, 0, 0, 2*pi, 0};
    double mass = 5.972e24; // Earth mass.
    Solarsystem q;
    q.add_planet(mass, U0);

    // arma::vec u0 = q.get_U0();
    // u0.print();
    q.solve_system();
    // arma::vec initial = q.get_U0();
    return 0;
}
