#include "energy_solver.cpp"

// TEST_CASE("test_if_energy_for_four_given_configurations_match_analytic_answer_for_dim_2")
void test_if_calculated_energy_match_analytical_answer_for_dim_2()
{
    double total_energy;
    double total_magnetization;
    int n = 2;
    double spin1[4] = {1, 1, 1, 1};
    double spin2[4] = {1, -1, -1, 1};
    double spin3[4] = {1, 1, -1, 1};
    double spin4[4] = {-1, -1, -1, -1};
    double energy[4] = {-8, 8, 0, -8};

    // mat(dimension, configuration)
    CircularMatrix mat1(2, spin1);
    CircularMatrix mat2(2, spin2);
    CircularMatrix mat3(2, spin3);
    CircularMatrix mat4(2, spin4);

    REQUIRE(IsingModel total_energy_and_magnetization(spin1, n,
            total_energy, total_magnetization)) == -8;
    REQUIRE(IsingModel total_energy_and_magnetization(spin2, n,
            total_energy, total_magnetization)) == 8;
    REQUIRE(IsingModel total_energy_and_magnetization(spin3, n,
            total_energy, total_magnetization)) == 0;
    REQUIRE(IsingModel total_energy_and_magnetization(spin4, n,
            total_energy, total_magnetization)) == -8;

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