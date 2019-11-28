#ifndef VECTOR_3D_H
#define VECTOR_3D_H

#include <cmath>
#include <iostream>

class Vector_3D
{
public:
    double* vector = new double[3];
    double x;
    double y;
    double z;

    Vector_3D(double x_input, double y_input, double z_input);
    void print();
    double& operator[] (int xyz);
    double dot(Vector_3D vec2);
    double abs_val();
    Vector_3D operator+ (Vector_3D vec2);
    Vector_3D operator- (Vector_3D vec2);
    Vector_3D operator* (double scalar);
    ~Vector_3D();
};

#endif