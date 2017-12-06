#include <iostream>
#include <iomanip>

using namespace std;

#include "Vector.hpp"
#include "rayTracer.hpp"
#include "geometry.hpp"

Sphere light_source(Point(10,10,10),3);
Point camera(0,3,10);
Sphere sphere(Point(0,0,0),2);
Plane table(Point(0,-2,0),Point(0,1,0) - Point(0,0,0));

int trace(Line& ray, int remaining){
    Point p = sphere.intersect(ray);
    int r = 0x87;
    int g = 0xce;
    int b = 0xeb;
    double ref;
    if (p.is_valid()){
	ray.reflect(p, (p-sphere.center).normalize());
	r = g = b = 0xff;
	ref = 1;
    }
    else{
	p = light_source.intersect(ray);
	if(p.is_valid()){
	    ray.reflect(p, (p-sphere.center).normalize());
	    r = g = b = 0xff;
	    ref = 0;
	}
	else{
	    p = table.intersect(ray);
	    if(p.is_valid()){
		int c = table.color(p);
		r = (c>>16)&0xff;
		g = (c>>8)&0xff;
		b = c&0xff;
		ray.reflect(p, table.normal);
		ref = 0.5;
	    }
	}
    }
    if (p.is_valid()){
	int c = (r<<16)+(g<<8)+b;
	if (remaining !=0){
	    ray.direction = ray.direction*-1;
	    c = trace(ray, remaining-1);
	}
	r = (1-ref)*r+ref*((c>>16)&0xff);
	g = (1-ref)*g+ref*((c>>8)&0xff);
	b = (1-ref)*b+ref*(c&0xff);
	int shadows = 0;
	for(int i = 0; i<SHADOW_SAMPLES; i++){
	    double theta = (double)rand()/RAND_MAX*M_PI;
	    double phi = (double)rand()/RAND_MAX*2*M_PI;
	    double r = (double)rand()/RAND_MAX*light_source.radius;
	    double x = r*sin(theta)*cos(phi) + light_source.center.x;
	    double y = r*sin(theta)*sin(phi) + light_source.center.y;
	    double z = r*cos(theta) + light_source.center.z;
	    Line shadow_ray(p, Point(x,y,z));
	    if (sphere.intersect(shadow_ray).is_valid()) shadows++;
	}
	double shadow_percent = 1-(double)shadows/SHADOW_SAMPLES;
	r*=shadow_percent;
	g*=shadow_percent;
	b*=shadow_percent;
    }
    return (r<<16)+(g<<8)+b;
}

int main(int argc, char** argv){
    int img[WIDTH*HEIGHT];
    for (int y = 0; y < HEIGHT; y++){
	for (int x = 0; x < WIDTH; x++) {
	    int r=0, g=0, b=0;
	    for (int i = 0; i<PIXEL_SAMPLES; i++){
		double xShift = (double)rand()/RAND_MAX-0.5;
		double yShift = (double)rand()/RAND_MAX-0.5;
		Point pixel((double)(x-WIDTH/2+xShift)/SCALE,(double)(HEIGHT/2-y+yShift)/SCALE,0);
		Line ray = Line(camera, pixel);
		int c = trace(ray, DEPTH);
		r += (c>>16)&0xff;
		g += (c>>8)&0xff;
		b += c&0xff;
	    }
	    r/=PIXEL_SAMPLES;
	    g/=PIXEL_SAMPLES;
	    b/=PIXEL_SAMPLES;
	    img[y*WIDTH+x]=(r<<16)+(g<<8) +b;
	    cout << "\r" << "Progress: " << fixed << setprecision(2) << ((double)((y*WIDTH)+x)/WIDTH/HEIGHT)*100 << "%" << flush;
	}
    }
    cout << endl;
    writeImage((char*)"test.png", WIDTH, HEIGHT, img);
    return 0;
}
