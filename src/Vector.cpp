#include <iostream>
#include "Vector.hpp"

double& Vector::operator[] (int n){
    return nums[n];
}
double Vector::operator[] (int n) const{
    return nums[n];
}
Vector Vector::operator+(const Vector& other) const{
    Vector vec;
    for(int i = 0; i < 3; i++){
	vec[i] = nums[i] + other[i];
    }
    return vec;
}
Vector Vector::operator*(double s) const{
    Vector vec;
    for(int i = 0; i < 3; i++){
	vec[i] = s* nums[i];
    }
    return vec;
}
Vector Vector::operator-(const Vector& other) const{
    return operator+(other*-1);
}
Vector& Vector::operator=(const Vector& other){
    for (int i = 0; i<3; i++){
	nums[i] = other[i];
    }
    return *this;
}
double Vector::dot(const Vector& other) const{
    double total = 0;
    for (int i = 0; i < 3; i++){
	total += nums[i]*other[i];
    }
    return total;
}
double Vector::length() const{
    return sqrt(dot(*this));
}
Vector Vector::cross(const Vector& other) const{
    Vector vec;
    vec[0] = nums[1]*other[2] - nums[2]*other[1];
    vec[1] = nums[2]*other[0] - nums[0]*other[2];
    vec[2] = nums[0]*other[1] - nums[1]*other[0];
    return vec;
}
Vector& Vector::normalize(){
    double l = length();
    for (int i = 0; i< 3; i++){
	nums[i] /=l;
    }
    return *this;
}

std::ostream& operator<<(std::ostream& strm, const Vector& v) {
    strm << "<" <<  v[0];
    for(int i = 1; i<3; i++){
	strm << "," << v[i];
    }
    strm << ">";
    return strm;
}

Vector operator*(double s, const Vector& v){
    return v*s;
}
