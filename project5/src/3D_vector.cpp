#include "vector_3D.h"

Vector_3D::Vector_3D(double x_input, double y_input, double z_input)
{
    vector[0] = x_input;
    vector[1] = y_input;
    vector[2] = z_input;

    x = x_input;
    y = y_input;
    z = z_input;
}

void Vector_3D::print()
{   /*
    Prints a nice visualization of the 3D-vector.
    */
    std::cout << "[" << vector[0] << vector[1] << vector[2] << "]" << std::endl;
    // std::cout << "[" << x << y << z << "]" << std::endl;
}

double& Vector_3D::operator[] (int xyz)
{   /*
    Index the matrix. Translate 2D indices to flat indices.

    Parameters
    ----------
    xyz : int
        Index of 3D-vector.
    */
    if ((xyz > 2) or (xyz < 0))
    {
        std::cout << "Index Error: Index out of bounds. Allowed indecies are 0, 1 and 2, not " << xyz << std::endl; 
        exit(EXIT_FAILURE);
    }
    // if (xyz == 0) {return x;}
    // else if (xyz == 1) {return y;}
    // else if (xyz == 2) {return z;}
    return vector[xyz];
}

double Vector_3D::dot(Vector_3D vec2)
{
    return vector[0]*vec2[0] + vector[1]*vec2[1] + vector[2]*vec2[2];
    // return x*vec2[0] + y*vec2[1] + z*vec2[2];
}

double abs_val()
{
    return std::sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
    // return std::sqrt(this.dot(this));
    // return std::sqrt(x*x + y*y + z*z);
}

Vector_3D Vector_3D::operator+ (Vector_3D vec2)
{
    return Vector_3D(vector[0]+vec2[0], vector[1]+vec2[1], vector[2]+vec2[2]);
}

Vector_3D Vector_3D::operator- (Vector_3D vec2)
{
    return Vector_3D(vector[0]-vec2[0], vector[1]-vec2[1], vector[2]-vec2[2]);
}

Vector_3D Vector_3D::operator* (double scalar)
{
    return Vector_3D(vector[0]*scalar, vector[1]*scalar, vector[2]*scalar);
}

Vector_3D::~Vector_3D()
{
    delete[] vector;
}