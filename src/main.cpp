#include <iostream>
#include <iomanip>
#include <thread>
#include <mutex>
#include <chrono>
#include <random>

using namespace std;

#include "Vector.hpp"
#include "rayTracer.hpp"
#include "geometry.hpp"

Sphere light_source(Point(20,20,20),3);
Point camera(0,3,10);
Sphere sphere1(Point(-3,-1,0),1);
Sphere sphere2(Point(0,-1,0),1);
Sphere sphere3(Point(3,-1,0),1);
Sphere sphere4(Point(0,2,-6),4);
Sphere sphere5(Point(0,13,-30),15);
Plane table(Point(0,-2,0),Point(0,1,0) - Point(0,0,0));
#define NUM_OBJS 7
Solid* objs[NUM_OBJS] {&light_source, &sphere1, &sphere2, &sphere3, &sphere4, &sphere5, &table};
default_random_engine generators[NUM_THREADS];
long num_rays[NUM_THREADS];

int trace(const Line& ray, int remaining, int thread_num){
    num_rays[thread_num]++;
    Point p(INFINITY,INFINITY,INFINITY);
    int r = 0x87;
    int g = 0xce;
    int b = 0xeb;
    Material mat = {0,0,0,0,0};
    Vector<3> normal;
    for(int i = 0; i < NUM_OBJS; i++){
	Point p0 = objs[i]->intersect(ray);
	if (p0.is_valid() && ((p-ray.point).length() > (p0-ray.point).length())){
	    p = p0;
	    normal = objs[i]->normal(p);
	    mat = objs[i]->material;
	    int c = objs[i]->get_color(p);
	    r = (c>>16)&0xff;
	    g = (c>>8)&0xff;
	    b = c&0xff;
	}
    }
    if (p.is_valid()){
	int c = (r<<16)+(g<<8)+b;
	if ((remaining!=0) && (!mat.is_light) && (mat.ref!=0)){
	    int refR = 0, refG=0, refB=0;
	    Line new_ray = ray;
	    new_ray.reflect(p, normal);
	    new_ray.direction = new_ray.direction*-1;
	    Vector<3> copy = new_ray.direction;
	    const int num = (mat.scatter_angle!=0)?SCATTER_SAMPLES:1; 
	    for(int i = 0; i< num; i++){
		Vector<3> v;
		double theta = (double)generators[thread_num]()/generators[thread_num].max()*mat.scatter_angle;
		double phi = (double)generators[thread_num]()/generators[thread_num].max()*2*M_PI;
		v[0] = sin(theta)*cos(phi);
		v[1] = sin(theta)*sin(phi);
		v[2] = cos(theta);
		if (fabs(copy[2] - 1) < 0.0001){
		    new_ray.direction = v;
		}
		else{
		    theta = acos(copy[2]);
		    Vector<3> n;
		    n[0]=0;n[1]=0;n[2]=1;
		    n=n.cross(copy);
		    n.normalize();
		    double sint = sin(theta);
		    double cost = cos(theta);
		    new_ray.direction = cost*(v-n*(n.dot(v))) + n*(n.dot(v)) + sint*n.cross(v);
		}
		c = trace(new_ray, remaining-1, thread_num);
		refR += (c>>16)&0xff;
		refG += (c>>8)&0xff;
		refB += (c)&0xff;
	    }
	    refR/=num;
	    refG/=num;
	    refB/=num;
	    c = (refR<<16)+(refG<<8)+refB;
	}
	int shadows = 0;
	for(int i = 0; i<SHADOW_SAMPLES; i++){
	    double theta = (double)generators[thread_num]()/generators[thread_num].max()*M_PI;
	    double phi = (double)generators[thread_num]()/generators[thread_num].max()*2*M_PI;
	    double r = (double)generators[thread_num]()/generators[thread_num].max()*light_source.radius;
	    double x = r*sin(theta)*cos(phi) + light_source.center.x;
	    double y = r*sin(theta)*sin(phi) + light_source.center.y;
	    double z = r*cos(theta) + light_source.center.z;
	    Line shadow_ray(p, Point(x,y,z));
	    for(int j = 1; j < NUM_OBJS; j++){
		if (objs[j]->intersect(shadow_ray).is_valid()){
		    shadows++;
		    break;
		}
	    }
	}
	double shadow_percent = (1-(double)shadows/SHADOW_SAMPLES);
	Line ref (light_source.center,p);
	ref.reflect(p,normal);
	double specular = pow(ray.direction.dot(ref.direction),1);
	double diffuse = -ref.direction.dot(normal);
	if (specular<0) specular = 0;
	if (diffuse<0) diffuse = 0;
	if (!mat.is_light){
	    r = r*((1-mat.ref)*((specular*mat.specular+diffuse*mat.diffuse)*shadow_percent*((light_source.color>>16)&0xff)/0xff/*+0.529411*/) + mat.ref*((c>>16)&0xff)/0xff);
	    g = g*((1-mat.ref)*((specular*mat.specular+diffuse*mat.diffuse)*shadow_percent*((light_source.color>>8)&0xff)/0xff/*+0.807843*/) + mat.ref*((c>>8)&0xff)/0xff);
	    b = b*((1-mat.ref)*((specular*mat.specular+diffuse*mat.diffuse)*shadow_percent*((light_source.color)&0xff)/0xff/*+0.921568*/) + mat.ref*((c)&0xff)/0xff);
	}
    }
    return (r<<16)+(g<<8)+b;
}

int completed = 0;
mutex completed_lock;
void do_rays_i(int* img, int num){
    int local_count = 0;
    for (int y = num; y < HEIGHT; y+=NUM_THREADS){
	for (int x = 0; x < WIDTH; x++) {
	    int r=0, g=0, b=0;
	    for (int i = 0; i<PIXEL_SAMPLES; i++){
		double xShift = (double)generators[num]()/generators[num].max()-0.5;
		double yShift = (double)generators[num]()/generators[num].max()-0.5;
		Point pixel((double)(x-WIDTH/2+xShift)/SCALE,(double)(HEIGHT/2-y+yShift)/SCALE,0);
		Line ray = Line(camera, pixel);
		int c = trace(ray, DEPTH, num);
		r += (c>>16)&0xff;
		g += (c>>8)&0xff;
		b += c&0xff;
	    }
	    r/=PIXEL_SAMPLES;
	    g/=PIXEL_SAMPLES;
	    b/=PIXEL_SAMPLES;
	    img[y*WIDTH+x]=(r<<16)+(g<<8) +b;
	}
	if (completed_lock.try_lock()){
	    completed += local_count +1;
	    completed_lock.unlock();
	    local_count = 0;
	}
	else{
	    local_count++;
	}
	if (num==0) cout << "\r" << "Progress: " << fixed << setprecision(2) << ((double)(completed)/HEIGHT)*100 << "%" << flush;
    }
}

int main(int argc, char** argv){
    cout << "Scale:          " << SCALE << endl;
    cout << "Width:          " << WIDTH << endl;
    cout << "Height:         " << HEIGHT << endl;
    cout << "Shadow samples: " << SHADOW_SAMPLES << endl;
    cout << "Pixel samples:  " << PIXEL_SAMPLES << endl;
    cout << "Scatter rays:   " << SCATTER_SAMPLES << endl;
    cout << "Threads:        " << NUM_THREADS << endl;
    cout << "Depth:          " << DEPTH << endl;
    cout << "Max rays:       " << WIDTH*HEIGHT*PIXEL_SAMPLES*DEPTH*pow(SCATTER_SAMPLES,DEPTH) << endl;
    light_source.material = {0,true,0,0,0};
    light_source.color = 0xffffff;
    sphere1.material = {0,false,0,0,1};
    sphere1.color = 0xffa0a0;
    sphere2.material = {1,false,0,1,0};
    sphere2.color = 0xa0ffa0;
    sphere3.material = {0.9,false,0.2,0.9,0.1};
    sphere3.color = 0xa0a0ff;
    sphere4.material = {0.5,false,0,0.5,0.5};
    sphere4.color = 0xffffff;
    sphere5.material = {0.5,false,0,0.5,0.5};
    sphere5.color = 0xffff7f;
    table.material = {0.2,false,0.2,0,1};
    int img[WIDTH*HEIGHT];
    thread threads[NUM_THREADS-1];
    auto start = chrono::system_clock::now();
    for(int i = 1; i<NUM_THREADS; i++){
	generators[i].seed(0);
	threads[i-1] = thread(do_rays_i,img,i);
	//threads[i-1].join();
    }
    do_rays_i(img,0);
    for(int i = 1; i<NUM_THREADS; i++){
	threads[i-1].join();
    }
    auto end = chrono::system_clock::now();
    cout << endl;
    chrono::duration<double> elapsed = end-start;
    double sec = elapsed.count();
    int min = sec/60;
    sec = sec - min*60;
    long rays = 0;
    for (int i = 0; i< NUM_THREADS; i++){
	rays+= num_rays[i];
    }
    cout << "Time: " << min << "m" << sec << "s" << endl;
    cout << "Rays: " << rays << endl;
    writeImage((char*)"test.png", WIDTH, HEIGHT, img);
    return 0;
}
