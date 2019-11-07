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
    q1.new_dim(4);
    q1.ordered_spin();
    q1.print();
    // std::cout << q1(0, 8, true) << std::endl;
}

// TEST_CASE("test_if_energy_for_four_given_configurations_match_analytic_answer_for_dim_2")
void test_if_energy_for_four_given_configurations_match_analytic_answer_for_dim_2()
{
    double init1[4] = {1, 1, 1, 1};
    double init2[4] = {1, -1, -1, 1};
    double init3[4] = {1, 1, -1, 1};
    double init4[4] = {-1, -1, -1, -1};

    // mat(dimension, configuration)
    CircularMatrix mat1(2, init1);
    CircularMatrix mat2(2, init2);
    CircularMatrix mat3(2, init3);
    CircularMatrix mat4(2, init4);

    // mat1.print();
    // std::cout << std::endl;
    // mat2.print();
    // std::cout << std::endl;
    // mat3.print();
    // std::cout << std::endl;
    // mat4.print();

    
}

// TEST_CASE("test_that_ordered_spin_produces_all_spin_up")
void test_that_ordered_spin_produces_all_spin_up()
{

}

// TEST_CASE("test_that_new_dim_produces_right_dimension")
void test_that_new_dim_produces_right_dimension()
{

}

// TEST_CASE("test_indexing")
void test_indexing()
{

}

// TEST_CASE("test_indexing_boundary_check")
void test_indexing_boundary_check()
{

}

int main()
{
    // test_print();
    test_if_energy_for_four_given_configurations_match_analytic_answer_for_dim_2();

    return 0;
}