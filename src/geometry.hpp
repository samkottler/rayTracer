#include <cmath>
#include "Vector.hpp"

class Point{
public:
    double x, y, z;
    Point(){}
    Point(double x1, double y1, double z1){
	x = x1;
	y = y1;
	z = z1;
    }
    bool is_valid(){
	if (isnan(x) || isnan(y) || isnan(z)) return false;
	return true;
    }
    Point operator+ (const Vector<3>& v) const{
	return Point(x+v[0], y+v[1], z + v[2]);
    }
    Vector<3> operator- (const Point& other) const{
	Vector<3> v;
	v[0] = x-other.x;
	v[1] = y-other.y;
	v[2] = z-other.z;
	return v;
    }
};

class Line{
public:
    Point point;
    Vector<3> direction;
    Line(const Point& p, const Vector<3>& d){
	point = p;
	direction = d;
	direction.normalize();
    }
    Line(const Point& p0, const Point& p1){
	point = p0;
	direction = p1-p0;
	direction.normalize();
    }
    Line& reflect(const Point& pNew, const Vector<3>& normal){
	point = pNew;
	Vector<3> w = normal.dot(direction)*normal;
	direction = 2*w - direction;
	return *this;
    }
};

class Plane{
private:
    Point point;
    Vector<3> normal;
public:
    Plane(const Point& p, const Vector<3>& n){
	point = p;
	normal = n;
	normal.normalize();
    }
    Point intersect(const Line& line) const{
	Vector<3> line_to_plane = point - line.point;
	double t = line_to_plane.dot(normal)/normal.dot(line.direction);
	//cout << t << endl;
	if (t<0) t = NAN;
	if (t == INFINITY || t == -INFINITY) t = NAN;
	return line.point + t*line.direction;
    }
    int color(const Point& p){
	Vector<3> v = p - point;
	int x = (int)(v[0])%2;
	int z = (int)(v[2])%2;
	if (v[0]<0) x=-x;
	if (v[2]<0) z=-z;
	int c = (x==z)?0:128;
	if (v[0]<0)
	    c = (x==z)?128:0;
	return (c<<16)+(c<<8) + c;
    }
};
    
class Sphere{
public:
    double radius;
    Point center;
    Sphere(const Point& p, double r){
	center = p;
	radius = r;
    };
    Point intersect(const Line& line) const{
	Vector<3> center_to_line = line.point - center;
	//cout << "(" << center_to_line[0] << "," << center_to_line[1] << "," << center_to_line[2] << ")" << endl;
 	double b = 2*line.direction.dot(center_to_line);
	double c = center_to_line.dot(center_to_line) - radius*radius;
	double s = sqrt(b*b - 4*c);
	double t = (-b - s)/2;
	if (t<0.000001) t = (-b+s)/2;
	if (t<0.000001) t = NAN;
	return line.point + t*line.direction;
    }
};
