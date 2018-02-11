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
    Color diffuse;
    Color specular;
    double specular_exp;
    Color refraction_diffuse;
    Color refraction_specular;
    double refraction_specular_exp;
    double refraction_index;
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
    Line& refract(const Point& pNew, const Vector& normal, double n1, double n2);
};

class Solid{
public:
    bool is_light;
    virtual Point intersect(const Line& line) const = 0;
    virtual const Material& get_material(const Point& p) const = 0;
    virtual Vector normal(const Point& p) const = 0;
};

class Plane: public Solid{
public:
    Point point;
    Vector norm;
    Material material1;
    Material material2;
    Plane(){};
    Plane(const Point& p, const Vector& n);
    Point intersect(const Line& line) const;
    const Material& get_material(const Point& p) const;
    Vector normal(const Point& p) const;
};
    
class Sphere: public Solid{
public:
    double radius;
    Point center;
    Material material;
    Color color;
    Sphere(const Point& p, double r);
    Point intersect(const Line& line) const;
    const Material& get_material(const Point& point) const;
    Vector normal(const Point& p) const;
};

class Face: public Plane{
public:
    Point* verts;
    int num;
    Material material;
    Face(const Vector& n, int num_verts, Point* v);
    Point intersect(const Line& line) const;
    const Material& get_material(const Point& p) const;
    Vector normal(const Point& p) const;
};
