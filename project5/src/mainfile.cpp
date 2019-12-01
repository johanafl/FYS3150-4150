#include "solarsystem.h"

int main()
{   
    arma::vec U1 = {1, 0, 0, 0, 2*pi, 0};
    arma::vec U2 = {0, 5.2, 0, 2.75, 0, 0};
    double mass1 = 5.972e24; // Earth mass.
    double mass2 = 1.898e27; // Earth mass.
    Solarsystem q;
    q.add_planet(mass1, U1);
    q.add_planet(mass2, U2);

    // arma::vec u0 = q.get_U0();
    // u0.print();
    q.solve_system();
    // arma::vec initial = q.get_U0();
    return 0;
}
