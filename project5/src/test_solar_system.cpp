#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "solar_system.h"

class TestClass : public SolarSystem
{
public:
    TestClass() : SolarSystem() {}
    arma::vec get_U0_test()
    {
        return get_U0();
    }
    void get_values(int& in_num_planets, int& in_array_size, double*& in_pos, double*& in_vel, double*& in_mass)
    {
        in_num_planets = num_planets;
        in_array_size = array_size;
        in_pos = pos;
        in_vel = vel;
        in_mass = mass;
    }
};

TEST_CASE("add_body")
{
    TestClass test;

    const double a_mass = 1e24;
    const double b_mass = 2e24;
    const double c_mass = 3e24;
    arma::vec a_initial = {1, 0, 0, 0, 2*pi, 0};
    arma::vec b_initial = {0, 1, 0, 0, 0, 2*pi};
    arma::vec c_initial = {0, 0, 1, 2*pi, 0, 0};

    int in_num_planets, in_array_size;
    double* in_pos;
    double* in_vel;
    double* in_mass;

    test.add_celestial_body(a_mass, a_initial);
    test.add_celestial_body(b_mass, b_initial);
    test.add_celestial_body(c_mass, c_initial);

    test.get_values(in_num_planets, in_array_size, in_pos, in_vel, in_mass);

    REQUIRE(in_num_planets == 3);

    REQUIRE(in_mass[0] == a_mass/solar_mass);
    REQUIRE(in_mass[1] == b_mass/solar_mass);
    REQUIRE(in_mass[2] == c_mass/solar_mass);

    for(int i=0; i<3; i++)
    {
        REQUIRE(in_pos[i] == a_initial(i));
        REQUIRE(in_pos[i + 3] == b_initial(i));
        REQUIRE(in_pos[i + 6] == c_initial(i));
        REQUIRE(in_vel[i] == a_initial(i + 3));
        REQUIRE(in_vel[i + 3] == b_initial(i + 3));
        REQUIRE(in_vel[i + 6] == c_initial(i + 3));
    }
}

TEST_CASE("get_u0")
{
    TestClass test;

    const double earth_mass = 5.972e24;
    arma::vec earth_initial = {1, 0, 0, 0, 2*pi, 0};
    test.add_celestial_body(earth_mass, earth_initial);

    arma::vec U0 = test.get_U0_test();

    for (int i=0; i<6; i++)
    {
        REQUIRE(earth_initial(i) == U0(i));
    }
}

TEST_CASE("resize")
{
    TestClass test;

    const double a_mass = 1e24;
    const double b_mass = 2e24;
    const double c_mass = 3e24;
    arma::vec a_initial = {1, 0, 0, 0, 2*pi, 0};
    arma::vec b_initial = {0, 1, 0, 0, 0, 2*pi};
    arma::vec c_initial = {0, 0, 1, 2*pi, 0, 0};

    int in_num_planets, in_array_size;
    double* in_pos;
    double* in_vel;
    double* in_mass;

    test.add_celestial_body(a_mass, a_initial);
    test.add_celestial_body(b_mass, b_initial);
    test.add_celestial_body(c_mass, c_initial);

    test.get_values(in_num_planets, in_array_size, in_pos, in_vel, in_mass);

    REQUIRE(in_array_size == 3*4);
}