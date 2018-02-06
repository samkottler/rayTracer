#include <cmath>
#include "Vector.hpp"

class Color{
public:
    double r,g,b;
    Color();
    Color(double r1, double g1, double b1);
    int to_int(double max);
    Color operator+(const Color& other) const;
    Color operator*(const Color& other) const;
    Color operator/(double num) const;
    Color operator*(double num) const;
    bool is_zero() const;
};

struct Material{
    Color ref;
    bool is_light;
    double scatter_angle;
    double diffuse;
    double specular;
    double specular_exp;
};

class Point{
public:
    double x, y, z;
    Point();
    Point(double x1, double y1, double z1);
    bool is_valid() const;
    Point operator+ (const Vector& v) const;
    Vector operator- (const Point& other) const;
};

class Line{
public:
    Point point;
    Vector direction;
    Line(const Point& p, const Vector& d);
    Line(const Point& p0, const Point& p1);
    Line& reflect(const Point& pNew, const Vector& normal);
};

class Solid{
public:
    Material material;
    virtual Point intersect(const Line& line) const = 0;
    virtual const Color& get_color(const Point& p) const = 0;
    virtual Vector normal(const Point& p) const = 0;
};

class Plane: public Solid{
public:
    Point point;
    Vector norm;
    Color color1;
    Color color2;
    Plane(const Point& p, const Vector& n);
    Point intersect(const Line& line) const;
    const Color& get_color(const Point& p) const;
    Vector normal(const Point& p) const;
};
    
class Sphere: public Solid{
public:
    double radius;
    Point center;
    Color color;
    Sphere(const Point& p, double r);
    Point intersect(const Line& line) const;
    const Color& get_color(const Point& point) const;
    Vector normal(const Point& p) const;
};
