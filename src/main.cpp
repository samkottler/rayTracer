#include <iostream>
#include <iomanip>
#include <thread>
#include <mutex>
#include <chrono>
#include <random>
#include <vector>
#include <string>

using namespace std;

#include "Vector.hpp"
#include "rayTracer.hpp"
#include "geometry.hpp"

vector<Solid*>* read_json_scene(string filename);

Point camera(0,3,10);

Color ambient;
vector<Solid*> objs;
int num_objs = 0;
int num_lights = 0;
default_random_engine generators[NUM_THREADS];
long num_rays[NUM_THREADS];
double exposure;
int pixel_samples;
int width,height,scale;

struct Trace_return{
    Color color;
    Point point;
};

void get_intersection(const Line& ray, Color* c, Material* mat, Vector<3>* normal, Point* p){
    for(int i = 0; i < num_objs; i++){
	Point p0 = objs[i]->intersect(ray);
	if (p0.is_valid() && ((*p-ray.point).length() > (p0-ray.point).length())){
	    *p = p0;
	    *normal = objs[i]->normal(p0);
	    *mat = objs[i]->material;
	    *c = objs[i]->get_color(p0);
	}
    }
}

Color get_direct_diffuse(const Line& ray, const Point& p, Vector<3> normal, int thread_num){
    Color total_diffuse;
    for(int k = 0; k< num_objs; k++){
	if (!objs[k]->material.is_light) continue;
	Sphere light_source = *((Sphere*)objs[k]);
	int shadows = 0;
	for(int i = 0; i<SHADOW_SAMPLES; i++){
	    double theta = (double)generators[thread_num]()/generators[thread_num].max()*M_PI;
	    double phi = (double)generators[thread_num]()/generators[thread_num].max()*2*M_PI;
	    double r = (double)generators[thread_num]()/generators[thread_num].max()*light_source.radius;
	    double x = r*sin(theta)*cos(phi) + light_source.center.x;
	    double y = r*sin(theta)*sin(phi) + light_source.center.y;
	    double z = r*cos(theta) + light_source.center.z;
	    Line shadow_ray(p, Point(x,y,z));
	    for(int j = 0; j < num_objs; j++){
		if (!objs[j]->material.is_light){
		    if (objs[j]->intersect(shadow_ray).is_valid()){
			shadows++;
			break;
		    }
		}
	    }

	}
	double shadow_percent = (1-(double)shadows/SHADOW_SAMPLES);
	Line ref (light_source.center,p);
	ref.reflect(p,normal);
	double specular = -pow(ray.direction.dot(ref.direction),1);
	double diffuse = -ref.direction.dot(normal);
	if (specular<0) specular = 0;
	if (diffuse<0) diffuse = 0;
	double dist = (p-light_source.center).length()/10;
	if (dist<1) dist = 1;
	total_diffuse = total_diffuse + light_source.color/dist/dist*diffuse*shadow_percent;
    }
    return total_diffuse;
}

Trace_return trace(const Line& ray, int remaining, int thread_num){
    num_rays[thread_num]++;
    Point p(INFINITY,INFINITY,INFINITY);
    Color c = ambient;
    Material mat = {Color(),0,0,0,0,0};
    Vector<3> normal;
    get_intersection(ray,&c,&mat,&normal,&p);
    //if (mat.is_light) return {c,p};
    if (p.is_valid() && (!mat.is_light)){
	Color ref_color(0,0,0);
	if ((remaining!=0) && (!mat.is_light) && !mat.ref.is_zero()){
	    Line new_ray = ray;
	    new_ray.reflect(p, normal);
	    new_ray.direction = new_ray.direction*-1;
	    //Vector<3> copy = new_ray.direction;
	    const int num = (mat.scatter_angle!=0)?SCATTER_SAMPLES:1;
	    double theta = acos(normal[2]);
	    Vector<3> n;
	    n[0]=0;n[1]=0;n[2]=1;
	    n=n.cross(normal);
	    n.normalize();
	    double sint = sin(theta);
	    double cost = cos(theta);
	    /*Line ref = new_ray;
	    ref.reflect(p,normal);
	    ref.direction=ref.direction*-1;*/
	    for(int i = 0; i< num; i++){
		Vector<3> v;
		double theta_rand = (double)generators[thread_num]()/generators[thread_num].max()*M_PI/2;
		double phi_rand = (double)generators[thread_num]()/generators[thread_num].max()*2*M_PI;
		v[0] = sin(theta_rand)*cos(phi_rand);
		v[1] = sin(theta_rand)*sin(phi_rand);
		v[2] = cos(theta_rand);
		if (fabs(normal[2] - 1) < 0.0001){
		    new_ray.direction = v;
		}
		else{
		    new_ray.direction = cost*(v-n*(n.dot(v))) + n*(n.dot(v)) + sint*n.cross(v);
		}
		//Trace_return deaper = trace(new_ray, remaining-1, thread_num);
		double specular = 0;
		Line ref = new_ray;
		ref.reflect(p,normal);
		ref.direction=ref.direction*-1;
		specular = pow(ref.direction.dot(ray.direction),mat.specular_exp);
		double diffuse = new_ray.direction.dot(normal);
		//cout<<specular<<endl;
		if (specular<0) specular = 0;
		if (diffuse<0) diffuse = 0;
		//double dist = (p-deaper.point).length()/10;
		//if (dist<1) dist = 1;
		//if (!deaper.point.is_valid()) dist = 1;
		double intensity = mat.diffuse*diffuse+mat.specular*specular;
		if (intensity<0.1){i--; continue;}
		Trace_return deaper = trace(new_ray, remaining-1, thread_num);
		ref_color = ref_color + deaper.color*intensity;
	    }
	    ref_color = ref_color/num;
	}
	c = c*ref_color;// + c*(ambient/8);
    }
    return {c,p};
}

int completed = 0;
mutex completed_lock;
void do_rays_i(int* img, int num){
    int local_count = 0;
    for (int y = num; y < height; y+=NUM_THREADS){
	for (int x = 0; x < width; x++) {
	    //int r=0, g=0, b=0;
	    Color c;
	    for (int i = 0; i<pixel_samples; i++){
		double xShift = (double)generators[num]()/generators[num].max()-0.5;
		double yShift = (double)generators[num]()/generators[num].max()-0.5;
		Point pixel((double)(x-width/2+xShift)/scale,(double)(height/2-y+yShift)/scale,0);
		Line ray = Line(camera, pixel);
		c = c + trace(ray, DEPTH, num).color;
	    }
	    c = c/pixel_samples;
	    if (x==0 && y==0){
		//cout<<c.r<<" "<<c.g<<" " <<c.b<<" " <<hex<< c.to_int(1.0/exposure)<<endl;
	    }
	    img[y*width+x]=c.to_int(1.0/exposure);
	}
	if (completed_lock.try_lock()){
	    completed += local_count +1;
	    completed_lock.unlock();
	    local_count = 0;
	}
	else{
	    local_count++;
	}
	if (num==0) cout << "\r" << "Progress: " << fixed << setprecision(2) << ((double)(completed)/height)*100 << "%" << flush;
    }
}

int main(int argc, char** argv){
    objs=*read_json_scene("scene.json");
    cout << "Scale:          " << scale << endl;
    cout << "Width:          " << width << endl;
    cout << "Height:         " << height << endl;
    cout << "Shadow samples: " << SHADOW_SAMPLES << endl;
    cout << "Pixel samples:  " << pixel_samples << endl;
    cout << "Scatter rays:   " << SCATTER_SAMPLES << endl;
    cout << "Threads:        " << NUM_THREADS << endl;
    cout << "Depth:          " << DEPTH << endl;
    cout << "Max rays:       " << width*height*pixel_samples*DEPTH*pow(SCATTER_SAMPLES,DEPTH) << endl;
    //objs=*read_json_scene("scene.json");
    int* img = new int[width*height];
    thread threads[NUM_THREADS-1];
    auto start = chrono::system_clock::now();
    for(int i = 1; i<NUM_THREADS; i++){
	generators[i].seed(0);
	threads[i-1] = thread(do_rays_i,img,i);
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
    writeImage((char*)"test.png", width, height, img);
    return 0;
}
