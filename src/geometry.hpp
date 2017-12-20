#include <cmath>
#include "Vector.hpp"

class Color{
public:
    double r,g,b;
    Color(){
	r = g = b = 0;
    }
    Color(double r1, double g1, double b1){
	r = r1;
	g = g1;
	b = b1;
    }
    int to_int(double max){
	//double l = log(max+1);
	int r_int = (1-(1.0/(r/max+1)))*255;
	int g_int = (1-(1.0/(g/max+1)))*255;
	int b_int = (1-(1.0/(b/max+1)))*255;
	if (r_int>255) r_int = 255;
	if (g_int>255) g_int = 255;
	if (b_int>255) b_int = 255;
	return (r_int<<16) + (g_int<<8) + b_int;
    }
    Color operator+(const Color& other) const{
	Color c(r+other.r, g+other.g, b+other.b);
	return c;
    }
    Color operator*(const Color& other) const{
	Color c(r*other.r, g*other.g, b*other.b);
	return c;
    }
    Color operator/(double num) const{
	Color c(r/num, g/num, b/num);
	return c;
    }
    Color operator*(double num) const{
	Color c(r*num, g*num, b*num);
	return c;
    }
    bool is_zero() const{
	if (r < 0.0001 && g < 0.0001 && b < 0.0001) return true;
	return false;
    }
};

struct Material{
    Color ref;
    bool is_light;
    double scatter_angle;
};

class Point{
public:
    double x, y, z;
    Point(){}
    Point(double x1, double y1, double z1){
	x = x1;
	y = y1;
	z = z1;
    }
    bool is_valid() const{
	if (isnan(x) || isnan(y) || isnan(z) || x==INFINITY || y==INFINITY || z==INFINITY) return false;
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

class Solid{
public:
    Material material;
    virtual Point intersect(const Line& line) const = 0;
    virtual Color get_color(const Point& p) const = 0;
    virtual Vector<3> normal(const Point& p) const = 0;
};

class Plane: public Solid{
public:
    Point point;
    Vector<3> norm;
    Plane(const Point& p, const Vector<3>& n){
	point = p;
	norm = n;
	norm.normalize();
    }
    Point intersect(const Line& line) const{
	Vector<3> line_to_plane = point - line.point;
	double t = line_to_plane.dot(norm)/norm.dot(line.direction);
	if (t<0.001) t = NAN;
	if (t == INFINITY || t == -INFINITY) t = NAN;
	return line.point + t*line.direction;
    }
    Color get_color(const Point& p) const{
	//return Color(0.5,0.9,0.9);
	Vector<3> v = p - point;
	int x = (int)(v[0])%2;
	int z = (int)(v[2])%2;
	if (v[0]<0) x=-x;
	if (v[2]<0) z=-z;
	double c = (x==z)?0.5:1;
	if ((v[0]<0) != (v[2]<0))
	    c = (x==z)?1:0.5;
	return Color(0.2,c,c);
    }
    Vector<3> normal(const Point& p) const{
	return norm;
    }
};
    
class Sphere: public Solid{
public:
    double radius;
    Point center;
    Color color;
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
	if (t<0.001) t = (-b+s)/2;
	if (t<0.001) t = NAN;
	//cout << t << endl;
	return line.point + t*line.direction;
    }
    Color get_color(const Point& point) const{
	return color;
    }
    Vector<3> normal(const Point& p) const{
	return (p-center).normalize();
    }
};
