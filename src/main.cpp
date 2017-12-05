#include <iostream>

using namespace std;

#include "Vector.hpp"
#include "rayTracer.hpp"
#include "geometry.hpp"

Point light_source(10,10,10);

int main(int argc, char** argv){
    int img[256*256];
    Point camera(0,2,10);
    Sphere s(Point(0,0,0),2);
    Plane table(Point(0,-2,0),Point(0,1,0) - Point(0,0,0));
    Vector<3> camera_to_pixel;
    for (int y = 0; y < 256; y++){
	for (int x = 0; x < 256; x++) {
	    Point pixel((double)(x-128)/SCALE,(double)(128-y)/SCALE,0);
	    Line ray = Line(camera, pixel);
	    Point p = s.intersect(ray);
	    if (p.is_valid()){
		//cout << p - Point(0,0,0) << endl;
		ray.reflect(p, (p-s.center).normalize());
		double brightness = ray.direction.dot((p - light_source).normalize());
		Line shadow_ray(p, light_source);
		int not_in_shadow = 1;
		if (s.intersect(shadow_ray).is_valid()) not_in_shadow = 0;
		if (brightness < 0) brightness = 0;
		int c = 256*brightness*not_in_shadow;
		img[y*256+x] = (c<<16)+(c<<8)+c;
	    }
	    else{
		p = table.intersect(ray);
		if(p.is_valid()){
		    Line shadow_ray(p, light_source);
		    int not_in_shadow = 1;
		    if (s.intersect(shadow_ray).is_valid()) not_in_shadow = 0;
		    img[y*256+x] = table.color(p)*not_in_shadow;
		}
	        else img[y*256+x] = 0xffffff;
	    }
	}
    }
    writeImage((char*)"test.png", 256, 256, img);
    return 0;
}
