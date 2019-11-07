#define CATCH_CONFIG_MAIN
#include "energy_solver.h"
#include "catch.hpp"

TEST_CASE("test_if_calculated_energy_match_analytic_answer_for_dim_2")
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
    CircularMatrix mat1(2, spin1);
    CircularMatrix mat2(2, spin2);
    CircularMatrix mat3(2, spin3);
    CircularMatrix mat4(2, spin4);
    CircularMatrix mat5(2, spin5);
    CircularMatrix mat6(2, spin6);

    IsingModel q(n, 0, 1337);

    q.total_energy_and_magnetization(mat1, n, total_energy,
                                              total_magnetization);
    REQUIRE(total_energy == -8);

    total_energy = 0;
    
    q.total_energy_and_magnetization(mat2, n, total_energy,
                                              total_magnetization);
    REQUIRE(total_energy == 0);
    
    total_energy = 0;

    q.total_energy_and_magnetization(mat3, n, total_energy,
                                              total_magnetization);
    REQUIRE(total_energy == 0);

    total_energy = 0;

    q.total_energy_and_magnetization(mat4, n, total_energy,
                                              total_magnetization);
    REQUIRE(total_energy == 8);

    total_energy = 0;

    q.total_energy_and_magnetization(mat5, n, total_energy,
                                              total_magnetization);
    REQUIRE(total_energy == 0);

    total_energy = 0;

    q.total_energy_and_magnetization(mat6, n, total_energy,
                                              total_magnetization);
    REQUIRE(total_energy == -8);
}

TEST_CASE("test_if_calculated_magnetization_match_analytic_answer_for_dim_2")
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
    CircularMatrix mat1(2, spin1);
    CircularMatrix mat2(2, spin2);
    CircularMatrix mat3(2, spin3);
    CircularMatrix mat4(2, spin4);
    CircularMatrix mat5(2, spin5);
    CircularMatrix mat6(2, spin6);

    IsingModel q(n, 0, 1337);

    q.total_energy_and_magnetization(mat1, n, total_energy,
                                              total_magnetization);
    REQUIRE(total_magnetization == 4);

    total_magnetization = 0;
    
    q.total_energy_and_magnetization(mat2, n, total_energy,
                                              total_magnetization);
    REQUIRE(total_magnetization == 2);
    
    total_magnetization = 0;

    q.total_energy_and_magnetization(mat3, n, total_energy,
                                              total_magnetization);
    REQUIRE(total_magnetization == 0);

    total_magnetization = 0;

    q.total_energy_and_magnetization(mat4, n, total_energy,
                                              total_magnetization);
    REQUIRE(total_magnetization == 0);

    total_magnetization = 0;

    q.total_energy_and_magnetization(mat5, n, total_energy,
                                              total_magnetization);
    REQUIRE(total_magnetization == -2);

    total_magnetization = 0;

    q.total_energy_and_magnetization(mat6, n, total_energy,
                                              total_magnetization);
    REQUIRE(total_magnetization == -4);
}

// int main()
// {
//     // test_print();
//     test_if_energy_for_four_given_configurations_match_analytic_answer_for_dim_2();

//     return 0;
// }