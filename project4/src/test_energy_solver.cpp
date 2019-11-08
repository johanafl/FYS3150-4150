#define CATCH_CONFIG_MAIN
#include "energy_solver.h"
#include "catch.hpp"


TEST_CASE("test_if_calculated_energy_and_magnetization_match_analytic_answer")
{
    double total_energy;
    double total_magnetization;
    int n = 2;
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

TEST_CASE("test_if_set_order_spins_gives_all_spins_up")
{
    int n = 7;

    IsingModel q(n, 0, 1337);

    q.set_order_spins();

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            REQUIRE(q.spin(i, j) == 1);
        }
    }
}

TEST_CASE("test_if_metropolis_flap_flips_a_value_or_not_depending_on_the_random_value")
{
    int n   = 2;
    int row = 0;
    int col = 1;
    double r1 = 1e4;
    double r2 = -1e4;
    double temp = 1;
    double total_energy;
    double total_magnetization;

    IsingModel q(n, 0, 1337);

    double* exp_delta_energy = new double[17];
    exp_delta_energy[0]  = std::exp(8*q.J/temp);
    exp_delta_energy[4]  = std::exp(4*q.J/temp);
    exp_delta_energy[8]  = 1;
    exp_delta_energy[12] = std::exp(-4*q.J/temp);
    exp_delta_energy[16] = std::exp(-8*q.J/temp);

    q.set_order_spins();
    q.metropolis_flap(q.spin, total_energy, total_magnetization, row, col, r1,
                      temp, exp_delta_energy);
    
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            REQUIRE(q.spin(i, j) == 1);
        }
    }

    q.set_order_spins();
    q.metropolis_flap(q.spin, total_energy, total_magnetization, row, col, r2,
                      temp, exp_delta_energy);
    
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            if (i == row and j == col)
            {
                REQUIRE(q.spin(i, j) == -1);
            }
            else
            {
                REQUIRE(q.spin(i, j) == 1);
            }
        }
    }
}