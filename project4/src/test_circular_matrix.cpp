#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "circular_matrix.h"


// TEST_CASE("test_print")
void test_print()
{   
    int seed = 1337;
    int n = 3;
    
    CircularMatrix q1(n);
    CircularMatrix q2(n, seed);
    q1.print();
    std::cout << std::endl;
    q2.print();
    std::cout << std::endl;
    q1.set_new_dim(4);
    q1.ordered_spin();
    q1.print();
    // std::cout << q1(0, 8, true) << std::endl;
}

TEST_CASE("test_that_ordered_spin_produces_all_spin_up")
{
    int n    = 2;
    int seed = 1337;

    CircularMatrix q(n, seed);

    q.ordered_spin();

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            REQUIRE(q(i, j) == 1);
    }
}

TEST_CASE("test_that_new_dim_and_seed_produces_new_dimension_and_seed")
{
    int n_old = 2;
    int n_new = 5;
    double old_seed = 1337;
    double new_seed = 2411;

    CircularMatrix mat(n_old, old_seed);

    mat.set_new_dim_and_seed(n_new, new_seed);

    REQUIRE(mat.dim == n_new);
    REQUIRE(mat.seed == new_seed);
}

TEST_CASE("test_indexing")
{
    int n    = 2;
    int seed = 1337;
    int row  = 0;
    int col  = 1;

    double spin1[4]  = {-1, -1, 1, -1};
    double spin2[4]  = {1, 1, -1, 1};

    CircularMatrix mat1(n, spin1);
    CircularMatrix mat2(n, spin2);

    REQUIRE(mat1(row, col, true));
    REQUIRE(mat2(row, col, true));
    REQUIRE(mat1(row, col));
    REQUIRE(mat2(row, col));
}

// TEST_CASE("test_indexing_boundary_check")
void test_indexing_boundary_check()
{
    // google how to require (test) EXIT_FAILURE.
}

// int main()
// {
//     // test_print();

//     return 0;
// }