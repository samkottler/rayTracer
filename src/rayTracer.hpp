#ifndef __MY_RAY_TRACER_HPP
#define __MY_RAY_TRACER_HPP

struct Trace_return{
    Color color;
    Point point;
};

extern "C" {
    int writeImage(char* filename, int width, int height, int* image);
}

vector<Solid*>* read_json_scene(string filename);
void get_intersection(const Line& ray, Color* c, Material* mat, Vector* normal, Point* p);
Color get_direct_radiance(const Line& ray, const Point& p, const Vector& normal, const Material& mat, int thread_num);
Trace_return trace(const Line& ray, int remaining, int thread_num);
void blur(Color* colors, int radius, double stddev);
void expose(int* img, Color* colors);
void do_rays_i(Color* img, int num);


#endif
