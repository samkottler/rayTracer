#include <iostream>
#include <iomanip>
#include <thread>
#include <mutex>
#include <chrono>
#include <random>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

#include "Vector.hpp"
#include "rayTracer.hpp"
#include "colors.hpp"
#include "geometry.hpp"

Sphere light_source(Point(-10,10,-10),3);
Point camera(0,3,10);
Sphere sphere1(Point(0,0,0),2);
Sphere sphere2(Point(0,-1,0),1);
Sphere sphere3(Point(3,-1,0),1);
Sphere sphere4(Point(0,2,-6),4);
Plane table(Point(0,-2,0),Point(0,1,0) - Point(0,0,0));
#define NUM_OBJS 3
Solid* objs[NUM_OBJS] {&light_source, &sphere1, &table};
default_random_engine generators[NUM_THREADS];

Spectrum trace(const Line& ray, int remaining, int thread_num){
    Point p(INFINITY,INFINITY,INFINITY);
    Spectrum to_return;
    for(int i = 0; i< 440/SPECTRUM_RESOLUTION; i++){
	//if (i<110/SPECTRUM_RESOLUTION)
	//to_return.samples[i] = 1-pow((double)i/440*SPECTRUM_RESOLUTION,0.5);
	//else
	to_return.samples[i] = 0.1;
    }
    Material mat = {0,0,0};
    Vector<3> normal;
    for(int i = 0; i < NUM_OBJS; i++){
	Point p0 = objs[i]->intersect(ray);
	if (p0.is_valid() && ((p-ray.point).length() > (p0-ray.point).length())){
	    p = p0;
	    normal = objs[i]->normal(p);
	    mat = objs[i]->material;
	    to_return = objs[i]->get_color(p);
	}
    }
    if (p.is_valid()){
	Spectrum reflected; 
	for (int i = 0; i< 440/SPECTRUM_RESOLUTION; i++) reflected.samples[i] = to_return.samples[i];
	if (remaining !=0){
	    for (int i = 0; i< 440/SPECTRUM_RESOLUTION; i++) reflected.samples[i] = 0;
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
		reflected += trace(new_ray, remaining-1, thread_num);
	    }
	    reflected/=num;
	    //reflected/=(DEPTH-remaining+1);
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
	if (mat.is_light) shadow_percent = 1;
	Spectrum light = light_source.color;
	light*=shadow_percent;
	if (mat.is_light) light/=(DEPTH-remaining+1);
	if ((remaining!=0)&!mat.is_light){
	  light+=reflected;
	}
	to_return*=light;
	//to_return/=remaining;
	//if (to_return.samples[0]>1) cout<< to_return.samples[0] << p-Point(0,0,0)<<endl;
    }
    return to_return;
}

int completed = 0;
mutex completed_lock;
void do_rays_i(XYZ* colors, int num){
    int local_count = 0;
    for (int y = num; y < HEIGHT; y+=NUM_THREADS){
	for (int x = 0; x < WIDTH; x++) {
	    Spectrum c;
	    for (int i = 0; i<PIXEL_SAMPLES; i++){
		double xShift = (double)generators[num]()/generators[num].max()-0.5;
		double yShift = (double)generators[num]()/generators[num].max()-0.5;
		Point pixel((double)(x-WIDTH/2+xShift)/SCALE,(double)(HEIGHT/2-y+yShift)/SCALE,0);
		Line ray = Line(camera, pixel);
		c += trace(ray, DEPTH, num);
	    }
	    c/=PIXEL_SAMPLES;
	    colors[y*WIDTH+x]=c.to_XYZ();
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
    light_source.material = {0,true,0};
    sphere1.material = {0.9,false,0};
    /*sphere1.color = 0xffa0a0;
    sphere2.material = {0.9,false,0};
    sphere2.color = 0xa0ffa0;
    sphere3.material = {0.9,false,0};
    sphere3.color = 0xa0a0ff;
    sphere4.material = {0.7,false,0.2};
    sphere4.color = 0xffffff;*/
    for(int i = 0; i<440/SPECTRUM_RESOLUTION; i++){
	light_source.color.samples[i] = 1;
	sphere1.color.samples[i] = 1;
    }
    table.material = {0.5,false,0};
    XYZ colors[WIDTH*HEIGHT];
    int img[WIDTH*HEIGHT];
    thread threads[NUM_THREADS-1];
    auto start = chrono::system_clock::now();
    for(int i = 1; i<NUM_THREADS; i++){
	generators[i].seed(0);
	threads[i-1] = thread(do_rays_i,colors,i);
    }
    do_rays_i(colors,0);
    for(int i = 1; i<NUM_THREADS; i++){
	threads[i-1].join();
    }
    double maxY = 0;
    int iSave = 0;
    for(int i = 0; i<HEIGHT*WIDTH; i++){
	if (maxY<colors[i].y){
	    maxY = colors[i].y;
	    //cout << i/WIDTH << "," << i%WIDTH<<endl;
	    iSave = i;
	}
    }
    for(int i = 0; i<HEIGHT*WIDTH; i++){
	colors[i].normalize(maxY);
	RGB color = colors[i].to_RGB();
	img[i] = color.to_int();
	//cout << colors[i].x << "," << colors[i].y << "," << colors[i].z << " " << color.r << "," << color.g << "," << color.b << " " << img[i] << endl;
    }
    img[iSave] =0xff0000;
    auto end = chrono::system_clock::now();
    cout << endl;
    chrono::duration<double> elapsed = end-start;
    double sec = elapsed.count();
    int min = sec/60;
    sec = sec - min*60;
    cout << "Time: " << min << "m" << sec << "s" << endl;
    writeImage((char*)"test.png", WIDTH, HEIGHT, img);
    return 0;
}
