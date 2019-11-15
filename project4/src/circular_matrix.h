#ifndef CIRCULAR_MATRIX_H
#define CIRCULAR_MATRIX_H

#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <cmath>
#include <string>


class CircularMatrix
{
public:
    int dim;
    double seed;
    double* matrix;

    CircularMatrix(int n, double seed_input);
    CircularMatrix(int n);
    CircularMatrix(int n, double* init_set);
    void initial_spin();
    void initial_spin(double seed_input);
    void ordered_spin();
    void set_new_dim_and_seed(int n, double new_seed);
    void set_new_dim(int n);
    void print();
    double& operator() (int row, int col);
    double& operator() (int row, int col, bool safe);
    ~CircularMatrix();
};

#endif