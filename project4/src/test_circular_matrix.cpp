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

    return 0;
}