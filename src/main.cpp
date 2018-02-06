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
#include "geometry.hpp"
#include "rayTracer.hpp"

Point camera(0,3,10);

Color ambient;
vector<Solid*> objs;
int num_objs = 0;
int num_lights = 0;
default_random_engine* generators;
long* num_rays;
double exposure;
int num_blurs;
int pixel_samples;
int width,height,scale;
double lens_radius;
double focal_length;
int shadow_samples;
int ray_depth;
int num_threads;
int scatter_samples;
double distance_scale;

void get_intersection(const Line& ray, Color* c, Material* mat, Vector* normal, Point* p){
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

Color get_direct_radiance(const Line& ray, const Point& p, Vector normal, Material mat, int thread_num){
    Color total_light;
    for(int k = 0; k< num_objs; k++){
	if (!objs[k]->material.is_light) continue;
	Sphere light_source = *((Sphere*)objs[k]);
	int shadows = 0;
	for(int i = 0; i<shadow_samples; i++){
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
	double shadow_percent = (1-(double)shadows/shadow_samples);
	Line ref (light_source.center,p);
	ref.reflect(p,normal);
	double specular = mat.specular*pow(ray.direction.dot(ref.direction),mat.specular_exp);
	double diffuse = -mat.diffuse*ref.direction.dot(normal);
	if (specular<0) specular = 0;
	if (diffuse<0) diffuse = 0;
	double dist = (p-light_source.center).length()/distance_scale;
	if (dist<1) dist = 1;
	total_light = total_light + light_source.color/dist/dist*(diffuse+specular)*shadow_percent;
    }
    return total_light;
}

Trace_return trace(const Line& ray, int remaining, int thread_num){
    num_rays[thread_num]++;
    Point p(INFINITY,INFINITY,INFINITY);
    Color c = ambient;
    Material mat = {Color(),0,0,0,0,0};
    Vector normal;
    get_intersection(ray,&c,&mat,&normal,&p);
    if (p.is_valid() && (!mat.is_light)){
	Color ref_color(0,0,0);
	if ((remaining!=0) && (!mat.is_light) && !mat.ref.is_zero()){
	    Line new_ray = ray;
	    new_ray.reflect(p, normal);
	    new_ray.direction = new_ray.direction*-1;
	    const int num = (mat.scatter_angle!=0)?scatter_samples:1;
	    double theta = acos(normal[2]);
	    Vector n;
	    n[0]=0;n[1]=0;n[2]=1;
	    n=n.cross(normal);
	    n.normalize();
	    double sint = sin(theta);
	    double cost = cos(theta);
	    Line ref = ray;
	    ref.reflect(p,normal);
	    ref.direction=ref.direction*-1;
	    double limit =0.1;// 0.5*(double)generators[thread_num]()/generators[thread_num].max();
	    for(int i = 0; i< num; i++){
		Vector v;
		double u_rand = (double)generators[thread_num]()/generators[thread_num].max();
		double phi_rand = (double)generators[thread_num]()/generators[thread_num].max()*2*M_PI;
		double sin_rand = sqrt(1-u_rand*u_rand);
		v[0] = sin_rand*cos(phi_rand);
		v[1] = sin_rand*sin(phi_rand);
		v[2] = u_rand;
		if (fabs(normal[2] - 1) < 0.0001){
		    new_ray.direction = v;
		}
		else{
		    new_ray.direction = cost*(v-n*(n.dot(v))) + n*(n.dot(v)) + sint*n.cross(v);
		}
		double specular = 0;
		specular = pow(ref.direction.dot(new_ray.direction),mat.specular_exp);
		double diffuse = new_ray.direction.dot(normal);
		if (specular<0) specular = 0;
		if (diffuse<0) diffuse = 0;
		double intensity = mat.diffuse*diffuse+mat.specular*specular;
		if (intensity<limit){i--; continue;}
		Trace_return deaper = trace(new_ray, remaining-1, thread_num);
		ref_color = ref_color + deaper.color*intensity;
	    }
	    ref_color = ref_color/num;
	}
	c = c*(ref_color+get_direct_radiance(ray,p,normal,mat,thread_num));
    }
    return {c,p};
}

void blur(Color* colors, int radius, double stddev){
    double kernel[(2*radius+1)*(2*radius+1)];
    double denom = 2*stddev*stddev;
    for (int y = -radius; y<=radius; y++){
	for (int x = -radius; x<=radius; x++){
	    kernel[(y+radius)*(2*radius+1)+x+radius] = exp(-(x*x+y*y)/denom)/denom/M_PI;
	}
    }
    Color* blurs[1+num_blurs];
    blurs[0] = new Color[width*height];
    for (int i = 0; i<width*height; i++){
	double r,g,b;
	r=g=b=0;
	if (colors[i].r*exposure > 1) r = exposure;
	if (colors[i].g*exposure > 1) g = exposure;
	if (colors[i].b*exposure > 1) b = exposure;
	blurs[0][i] = Color(r,g,b);
    }
    for (int i = 1; i< 1+num_blurs; i++){
	blurs[i] = new Color[width*height];
	for (int y = 0; y<height; y++){
	    for(int x = 0; x<width; x++){
		Color c;
		for (int ky = -radius; ky<=radius; ky++){
		    for (int kx = -radius; kx<=radius; kx++){
			if ((y+ky)>=0 && (y+ky)<height && (x+kx)>=0 && (x+kx)<width){
			    c = c+blurs[i-1][(y+ky)*width+(x+kx)]*kernel[(ky+radius)*(2*radius+1)+kx+radius];
			}
		    }
		}
		blurs[i][y*width+x] = c;
	    }
	}
    }
    for (int i = 0; i<width*height; i++){
	colors[i] = colors[i] + blurs[num_blurs][i];
    }
    for (int i = 1; i<num_blurs+1; i++){
	delete blurs[i];
    }
}

void expose(int* img, Color* colors){
    blur(colors, 20, 3);
    for (int y=0; y<height; y++){
	for (int x=0; x<width; x++){
	    img[y*width+x]=colors[y*width+x].to_int(1.0/exposure);
	}
    }
}

int completed = 0;
mutex completed_lock;
void do_rays_i(Color* img, int num){
    int local_count = 0;
    for (int y = num; y < height; y+=num_threads){
	for (int x = 0; x < width; x++) {
	    Color c;
	    for (int i = 0; i<pixel_samples; i++){
		double xShift = (double)generators[num]()/generators[num].max()-0.5;
		double yShift = (double)generators[num]()/generators[num].max()-0.5;
		Point pixel((double)(x-width/2+xShift)/scale,(double)(height/2-y+yShift)/scale,0);
		double theta = (double)generators[num]()/generators[num].max()*2*M_PI;
		double rad = (double)generators[num]()/generators[num].max()*lens_radius;
		Vector lens_shift;
		lens_shift[0] = rad*cos(theta); lens_shift[1] = rad*sin(theta); lens_shift[2] = 0;
		Vector D = pixel-camera;
		D.normalize();
		Point C = camera+focal_length*D;
		Line ray = Line(camera+lens_shift, C);
		c = c + trace(ray, ray_depth, num).color;
	    }
	    c = c/pixel_samples;
	    img[y*width+x]=c;
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
    cout << "Shadow samples: " << shadow_samples << endl;
    cout << "Pixel samples:  " << pixel_samples << endl;
    cout << "Scatter rays:   " << scatter_samples << endl;
    cout << "Threads:        " << num_threads << endl;
    cout << "Depth:          " << ray_depth << endl;
    cout << "Max rays:       " << (double)width*height*pixel_samples*ray_depth*pow(scatter_samples,ray_depth) << endl;
    num_rays = new long[num_threads];
    generators = new default_random_engine[num_threads];
    int* img = new int[width*height];
    Color* colors = new Color[width*height];
    thread threads[num_threads-1];
    auto start = chrono::system_clock::now();
    for(int i = 1; i<num_threads; i++){
	generators[i].seed(0);
	threads[i-1] = thread(do_rays_i,colors,i);
    }
    do_rays_i(colors,0);
    for(int i = 1; i<num_threads; i++){
	threads[i-1].join();
    }
    auto end = chrono::system_clock::now();
    cout << endl;
    chrono::duration<double> elapsed = end-start;
    double sec = elapsed.count();
    int min = sec/60;
    sec = sec - min*60;
    long rays = 0;
    for (int i = 0; i< num_threads; i++){
	rays+= num_rays[i];
    }
    cout << "Time: " << min << "m" << sec << "s" << endl;
    cout << "Rays: " <<defaultfloat<< (double)rays << endl;
    expose(img, colors);
    writeImage((char*)"test.png", width, height, img);
    return 0;
}
