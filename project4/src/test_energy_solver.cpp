#include "energy_solver.cpp"

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

int main()
{
    // test_print();
    test_if_energy_for_four_given_configurations_match_analytic_answer_for_dim_2();

    return 0;
}