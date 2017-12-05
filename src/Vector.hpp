#ifndef __MY_VECTOR_HPP
#define __MY_VECTOR_HPP
#include <cmath>
#include <iostream>

template <int N>
class Vector{
private:
    double nums[N];
public:
    const int size = N;
    double& operator[] (int n){
	return nums[n];
    }
    double operator[] (int n) const{
	return nums[n];
    }
    Vector<N> operator+ (const Vector<N>& other) const{
	Vector<N> vec;
	for(int i = 0; i < N; i++){
	    vec[i] = nums[i] + other[i];
	}
	return vec;
    }
    Vector<N> operator* (double s) const{
	Vector<N> vec;
	for(int i = 0; i < N; i++){
	    vec[i] = s* nums[i];
	}
	return vec;
    }
    Vector<N> operator- (const Vector<N>& other) const{
	return operator+(other*-1);
    }
    Vector<N>& operator= (const Vector<N>& other){
	for (int i = 0; i<N; i++){
	    nums[i] = other[i];
	}
	return *this;
    }
    double dot(const Vector<N>& other) const{
	double total = 0;
	for (int i = 0; i < N; i++){
	    total += nums[i]*other[i];
	}
	return total;
    }
    double length() const{
	return sqrt(dot(*this));
    }
    Vector<3> cross(const Vector<3>& other) const{
	Vector<3> vec;
	vec[0] = nums[1]*other[2] - nums[2]*other[1];
	vec[1] = nums[2]*other[0] - nums[0]*other[2];
	vec[2] = nums[0]*other[1] - nums[1]*other[0];
	return vec;
    }
    Vector<3>& normalize(){
	double l = length();
	for (int i = 0; i< N; i++){
	    nums[i] /=l;
	}
	return *this;
    }
};

template <int N>
std::ostream& operator<<(std::ostream& strm, const Vector<N>& v) {
    strm << "<" <<  v[0];
    for(int i = 1; i<v.size; i++){
	strm << "," << v[i];
    }
    strm << ">";
    return strm;
}

template <int N>
Vector<N> operator*(double s, const Vector<N>& v){
    return v*s;
}

#endif
