#define CATCH_CONFIG_MAIN
#include "energy_solver.h"
#include "catch.hpp"

void metropolis_flap(CircularMatrix& spin, double& total_energy,
        double& total_magnetization, int row, int col, double metropolis_random,
        double temperature, double* exp_delta_energy);

TEST_CASE("test_if_calculated_energy_and_magnetization_match_analytic_answer")
{
    double total_energy;
    double total_magnetization;
    int n = 2;
    double tol = 1e7;
    double spin1[4]  = {1, 1, 1, 1};
    double spin2[4]  = {1, 1, 1, -1};
    double spin3[4]  = {1, -1, 1, -1};
    double spin4[4]  = {1, -1, -1, 1};
    double spin5[4]  = {-1, -1, -1, 1};
    double spin6[4]  = {-1, -1, -1, -1};

    // mat(dimension, configuration)
    CircularMatrix mat1(n, spin1);
    CircularMatrix mat2(n, spin2);
    CircularMatrix mat3(n, spin3);
    CircularMatrix mat4(n, spin4);
    CircularMatrix mat5(n, spin5);
    CircularMatrix mat6(n, spin6);

    IsingModel q(n, 0, 1337);

    q.total_energy_and_magnetization(mat1, n, total_energy,
                                              total_magnetization);
    REQUIRE(total_energy == -8);
    REQUIRE(total_magnetization == 4);

    total_energy        = 0;
    total_magnetization = 0;
    q.total_energy_and_magnetization(mat2, n, total_energy,
                                              total_magnetization);
    REQUIRE(total_energy == 0);
    REQUIRE(total_magnetization == 2);
    
    total_energy        = 0;
    total_magnetization = 0;
    q.total_energy_and_magnetization(mat3, n, total_energy,
                                              total_magnetization);
    REQUIRE(total_energy == 0);
    REQUIRE(total_magnetization == 0);

    total_energy        = 0;
    total_magnetization = 0;
    q.total_energy_and_magnetization(mat4, n, total_energy,
                                              total_magnetization);
    REQUIRE(total_energy == 8);
    REQUIRE(total_magnetization == 0);

    total_energy        = 0;
    total_magnetization = 0;
    q.total_energy_and_magnetization(mat5, n, total_energy,
                                              total_magnetization);
    REQUIRE(total_energy == 0);
    REQUIRE(total_magnetization == -2);

    total_energy        = 0;
    total_magnetization = 0;
    q.total_energy_and_magnetization(mat6, n, total_energy,
                                              total_magnetization);
    REQUIRE(total_energy == -8);
    REQUIRE(total_magnetization == -4);
}

void set_new_input(int spin_mat_dim, int mc_iterations_input, double inter_strenght_J, long seed);

TEST_CASE("test_if_set_interactions_strength_gives_new_J_value")
{
    int n        = 2;
    double value = 4;

    IsingModel q(n, 0, 1337);

    q.set_interactions_strength(value);
    REQUIRE(q.J == value);
}

TEST_CASE("test_if_set_mc_iterations_gives_new_mc_iterations_value")
{
    int n        = 2;
    double value = 117;

    IsingModel q(n, 0, 1337);

    q.set_mc_iterations(value);
    REQUIRE(q.mc_iterations == value);
}

TEST_CASE("test_if_set_spin_dim_gives_new_dim_value")
{
    int n        = 2;
    double value = 17;

    IsingModel q(n, 0, 1337);

    q.set_spin_dim(value);
    REQUIRE(q.n == value);
}

TEST_CASE("test_if_set_order_spins_gives_new_dim_value")
{
    int n = 7;

    IsingModel q(n, 0, 1337);

    q.set_order_spins();

    for (int i = 0; i < 7; i++)
    {
        for (int j = 0; j < 7; j++)
        {
            REQUIRE(q.spin(i, j) == 1);
        }
    }
}

// int main()
// {
//     // test_print();
//     test_if_energy_for_four_given_configurations_match_analytic_answer_for_dim_2();

//     return 0;
// }