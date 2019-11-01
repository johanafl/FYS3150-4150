#ifndef CIRCULAR_MATRIX_H
#define CIRCULAR_MATRIX_H

#include <random>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <chrono>
#include <cmath>

const double pi = 3.14159265358979323846;

class CircularMatrix
{
private:
int dim;
int seed;

public:
    double* matrix;

    CircularMatrix(int n, double seed_input);
    CircularMatrix(int n);
    CircularMatrix(int n, double* init_set);
    void initial_spin();
    void ordered_spin();
    void print();
    double& operator() (int row, int col);
    double& operator() (int row, int col, bool safe);
    ~CircularMatrix();
};

#endif