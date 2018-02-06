#ifndef __MY_VECTOR_HPP
#define __MY_VECTOR_HPP
#include <cmath>
#include <iostream>

class Vector{
private:
    double nums[3];
public:
    double& operator[] (int n);
    double operator[] (int n) const;
    Vector operator+ (const Vector& other) const;
    Vector operator* (double s) const;
    Vector operator- (const Vector& other) const;
    Vector& operator= (const Vector& other);
    double dot(const Vector& other) const;
    double length() const;
    Vector cross(const Vector& other) const;
    Vector& normalize();
};

std::ostream& operator<<(std::ostream& strm, const Vector& v);
Vector operator*(double s, const Vector& v);

#endif
