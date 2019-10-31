#include "circular_matrix.h"

void test_print()
{   
    int seed = 1337;
    int n = 3;
    
    CircularMatrix q1(n);
    CircularMatrix q2(n, seed);
    q1.print();
    std::cout << std::endl;
    q2.print();
    
    std::cout << q1(0, 8, true) << std::endl;
}

void test_energy()
{
    double init1[4] = {1,1,1,1};
    double init2[4] = {1,-1,-1,1};
    double init3[4] = {1,1,-1,1};
    double init4[4] = {-1,-1,-1,-1};

    CircularMatrix mat1(2, init1);
    CircularMatrix mat2(2, init2);
    CircularMatrix mat3(2, init3);
    CircularMatrix mat4(2, init4);

    mat1.print();
    std::cout << std::endl;
    mat2.print();
    std::cout << std::endl;
    mat3.print();
    std::cout << std::endl;
    mat4.print();
}

int main()
{
    // test_print();
    test_energy();

    return 0;
}