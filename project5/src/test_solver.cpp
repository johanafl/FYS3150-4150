// #include "test_solver.h"
// #include <armadillo>
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "solver.h"

class TestClass
{
public:
    TestClass() {}
    arma::vec acceleration(arma::vec u, double t)
    {
        return arma::vec {0,0,0};
    }
};

TEST_CASE("set_initial_conditions")
{
    ForwardEuler<TestClass> q(1,1);
    arma::vec U0 = {1,2,3,4,5,6};
    q.set_initial_conditions(U0);
    arma::mat pos = q.get_pos();
    arma::mat vel = q.get_vel();
    for (int i=0; i<3; i++)
    {
        REQUIRE(pos(i,0) == U0(i));
        REQUIRE(vel(i,0) == U0(i + 3));
    }
}

TEST_CASE("euler")
{
    TestClass test;
    double dt = 1;
    arma::vec U0 = {1,0,0,1,0,0};
    
    ForwardEuler<TestClass> solved(10, 1);
    solved.set_initial_conditions(U0);
    solved.solve(test, dt);

    arma::mat pos = solved.get_pos();
    arma::mat vel = solved.get_vel();

    for (int i=0; i<10; i++)
    {
        REQUIRE(pos(0, i+1) == pos(0, i) + vel(0, i)*dt);
        REQUIRE(pos(1, i+1) == 0);
        REQUIRE(pos(2, i+1) == 0);
        REQUIRE(vel(0, i+1) == 1);
        REQUIRE(vel(1, i+1) == 0);
        REQUIRE(vel(2, i+1) == 0);
    }
}
