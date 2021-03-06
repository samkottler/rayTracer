#include <cmath>
#include <random>
#include "geometry.hpp"

using namespace std;

Color::Color(){
    r=b=g=0;
}
Color::Color(double r1, double g1, double b1){
    r = r1;
    g = g1;
    b = b1;
}
int Color::to_int(double max){
    int r_int = 255*r/max;
    int g_int = 255*g/max;
    int b_int = 255*b/max;
    if (r_int>255) r_int = 255;
    if (g_int>255) g_int = 255;
    if (b_int>255) b_int = 255;
    return (r_int<<16) + (g_int<<8) + b_int;
}
Color Color::operator+(const Color& other) const{
    Color c(r+other.r, g+other.g, b+other.b);
    return c;
}
Color Color::operator*(const Color& other) const{
    Color c(r*other.r, g*other.g, b*other.b);
    return c;
}
Color Color::operator/(double num) const{
    Color c(r/num, g/num, b/num);
    return c;
}
Color Color::operator*(double num) const{
    Color c(r*num, g*num, b*num);
    return c;
}
bool Color::is_zero() const{
    if (r < 0.0001 && g < 0.0001 && b < 0.0001) return true;
    return false;
}

Point::Point(){}
Point::Point(double x1, double y1, double z1){
    x = x1;
    y = y1;
    z = z1;
}
bool Point::is_valid() const{
    if (isnan(x) || isnan(y) || isnan(z) || x==INFINITY || y==INFINITY || z==INFINITY) return false;
    return true;
}
Point Point::operator+ (const Vector& v) const{
    return Point(x+v[0], y+v[1], z + v[2]);
}
Vector Point::operator- (const Point& other) const{
    Vector v;
    v[0] = x-other.x;
    v[1] = y-other.y;
    v[2] = z-other.z;
    return v;
}

Line::Line(){
}
Line::Line(const Point& p, const Vector& d){
    point = p;
    direction = d;
    direction.normalize();
}
Line::Line(const Point& p0, const Point& p1){
    point = p0;
    direction = p1-p0;
    direction.normalize();
}
Line& Line::reflect(const Point& pNew, const Vector& normal){
    point = pNew;
    Vector w = normal.dot(direction)*normal;
    direction = 2*w - direction;
    return *this;
}
Line& Line::refract(const Point& pNew, const Vector& normal, double n1, double n2){
    point = pNew;
    double ratio = n1/n2;
    double nv = fabs(normal.dot(direction));
    //cout<< nv << " " <<ratio*nv - sqrt(1-ratio*ratio*(1-nv*nv))<<endl;
    direction = ratio*direction + (ratio*nv - sqrt(1-ratio*ratio*(1-nv*nv)))*normal;
    direction.normalize();
    return *this;
}

Plane::Plane(const Point& p, const Vector& n){
    point = p;
    norm = n;
    norm.normalize();
}
Point Plane::intersect(const Line& line) const{
    Vector line_to_plane = point - line.point;
    double t = line_to_plane.dot(norm)/norm.dot(line.direction);
    if (t<0.001) t = NAN;
    if (t == INFINITY || t == -INFINITY) t = NAN;
    return line.point + t*line.direction;
}
const Material& Plane::get_material(const Point& p) const{
    Vector v = p - point;
    int x = (int)(v[0]/5)%2;
    int z = (int)(v[2]/5)%2;
    if (v[0]<0) x=-x;
    if (v[2]<0) z=-z;
    const Material* m = (x==z)?&material1:&material2;
    if ((v[0]<0) != (v[2]<0))
	m = (x==z)?(&material2):(&material1);
    return *m;
}
Vector Plane::normal(const Point& p) const{
    return norm;
}
Point Plane::rand_point(default_random_engine& gen) const {
    return point;
}


Sphere::Sphere(const Point& p, double r){
    center = p;
    radius = r;
};
Point Sphere::intersect(const Line& line) const{
    Vector center_to_line = line.point - center;
    double b = 2*line.direction.dot(center_to_line);
    double c = center_to_line.dot(center_to_line) - radius*radius;
    double s = sqrt(b*b - 4*c);
    double t = (-b - s)/2;
    if (t<0.001) t = (-b+s)/2;
    if (t<0.001) t = NAN;
    return line.point + t*line.direction;
}
const Material& Sphere::get_material(const Point& point) const{
    return material;
}
Vector Sphere::normal(const Point& p) const{
    return (p-center).normalize();
}
Point Sphere::rand_point(default_random_engine& gen) const{
    double u_rand = (double)gen()/gen.max();
    double phi_rand = (double)gen()/gen.max()*2*M_PI;
    double sin_rand = sqrt(1-u_rand*u_rand);
    return Point(center.x + sin_rand*cos(phi_rand), center.y + sin_rand*sin(phi_rand), center.z + u_rand);
}

Face::Face(const Vector& n, int num_verts, Point* v){
    point = v[0];
    norm = n;
    verts = v;
    num = num_verts;
    v1 = verts[0] - verts[1];
    v2 = verts[2] - verts[1];
}
Point Face::intersect(const Line& line) const{
    Point p = this->Plane::intersect(line);
    for (int i = 0; i<num; i++){
	double to_p = (verts[i]-verts[(i+1)%num]).dot(verts[i]-p);
	double to_next = (verts[i]-verts[(i+1)%num]).dot(verts[i]-verts[(i+2)%num]);
	if (to_p*to_next<0) return Point(NAN,NAN,NAN);
    }
    return p;
}
const Material& Face::get_material(const Point& p) const{
    return material;
}
Vector Face::normal(const Point& p) const{
    return norm;
}
Point Face::rand_point(default_random_engine& gen) const{
    double r1 = (double)gen()/gen.max();
    double r2 = (double)gen()/gen.max();
    return verts[1] + v1*r1 + v2*r2;
}

Path::Path(const Point& start){
    head = new path_node;
    head->point = start;
    head->next = nullptr;
}
Color Path::trace(){
    path_node* current = head->next;
    Color c = head->mat.ref;
    if (!head->point.is_valid()) return Color(0,0,0);
    Vector to_next = (head->point - current->point).normalize();
    while (current->next){
	Material mat = current->mat;
	Color diffuse;
	Color specular;
	double specular_exp;
	double s;
	if (current->reflection){
	    diffuse = mat.diffuse;
	    specular = mat.specular;
	    specular_exp = mat.specular_exp;
	    s = pow(current->comparison.dot(to_next),specular_exp);
	}
	else{
	    diffuse = mat.refraction_diffuse;
	    specular = mat.refraction_specular;
	    s = 1;
	}
	//double s = pow(current->comparison.dot(to_next),specular_exp);
	double d = current->normal.dot(to_next);
	if (s<0) s=0;
	if (d<0) d=0;
	c = c*(diffuse*d + specular*s);
	to_next = (current->point - current->next->point).normalize();
	current = current->next;
    }
    return c;
}
void Path::add(const Point p, const Vector normal, const Vector comparison, const Material mat, bool reflect){
    path_node* new_node = new path_node;
    new_node->next = head;
    new_node->point = p;
    new_node->normal = normal;
    new_node->comparison = comparison;
    new_node->mat = mat;
    new_node->reflection = reflect;
    head = new_node;
}
Path::~Path(){
    path_node* current = head;
    while (current){
	path_node* next = current->next;
	delete current;
	current = next;
    }
}
