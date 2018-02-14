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

Point camera(0,0,40);

Color ambient;
vector<Solid*> objs;
int num_objs = 0;
int num_lights = 0;
default_random_engine* generators;
long* num_rays;
double exposure;
int num_blurs;
int pixel_samples;
int width,height;
double scale;
double lens_radius;
double focal_length;
int shadow_samples;
int ray_depth;
int num_threads;
int scatter_samples;
double distance_scale;

void get_intersection(const Line& ray, Material* mat, Vector* normal, Point* p){
    for(int i = 0; i < num_objs; i++){
	Point p0 = objs[i]->intersect(ray);
	if (p0.is_valid() && ((*p-ray.point).length() > (p0-ray.point).length())){
	    *p = p0;
	    *normal = objs[i]->normal(p0);
	    *mat = objs[i]->get_material(p0);
	}
    }
}

Line get_indirect_reflection(const Line& ray, const Point& p, const Vector& normal, const Material& mat, int thread_num){
    Line new_ray(p,Vector());
    double theta = acos(normal[2]);
    Vector n;
    n[0]=0;n[1]=0;n[2]=1;
    n=n.cross(normal);
    n.normalize();
    double sint = sin(theta);
    double cost = cos(theta);
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
    return new_ray;
}			      

Line get_refraction(const Line& ray, const Point& p, const Vector& normal, const Material& mat, double refraction_index, int thread_num){
    double n2;
    Vector actual_normal = normal;
    if (fabs(mat.refraction_index-refraction_index) < 0.0001){ // leaving
	n2 = 1;
	actual_normal = -1*normal;
    }
    else // entering
	n2 = mat.refraction_index;
    Line ref = ray;
    double nv = actual_normal.dot(ray.direction);
    if((n2*n2)/(refraction_index*refraction_index) < 1-nv*nv) ref.reflect(p,actual_normal);
    else ref.refract(p,actual_normal,refraction_index,n2);
    return ref;
}

void add_to_path(const Line& ray, int remaining, double refraction_index, int thread_num, Path* path){
    //cout<<remaining<<endl;
    num_rays[thread_num]++;
    Point p(INFINITY,INFINITY,INFINITY);
    Material mat = {Color(),0,Color(),Color(),0};
    Vector normal;
    get_intersection(ray,&mat,&normal,&p);
    if (mat.is_light){
	path->add(p,Vector(), Vector(), mat, true);
	return;
    }
    else if (p.is_valid()){
	if (remaining == 0){
	    for(int i = 0; i<num_objs; i++){
		if (objs[i]->is_light){
		    Point rp = objs[i]->rand_point(generators[thread_num]);
		    for(int j = 0; j < num_objs; j++){
			if (i==j) continue;
			Point point = objs[j]->intersect(ray);
			if (point.is_valid() && (point-p).length()<(rp-p).length()){
			    path->add(Point(INFINITY,INFINITY,INFINITY),Vector(),Vector(),mat, true);
			    return;
			}
		    }
		    path->add(rp,Vector(),Vector(),objs[i]->get_material(rp),true);
		    return;
		}
	    }
	}
        Line forward;//=get_indirect_reflection(ray,p,normal,mat,thread_num);
	int change = 1;
	double n2 = 1;
	bool reflect = false;
	if (fabs(mat.refraction_index-refraction_index) < 0.0001){ // leaving
	    forward = get_refraction(ray,p,normal,mat,refraction_index,thread_num);
	    change = 0;
	}
	else{ // entering
	    if (mat.refraction_index == 0){
		forward = get_indirect_reflection(ray,p,normal,mat,thread_num);
		reflect = true;
	    }
	    else{
		forward = get_refraction(ray,p,normal,mat,refraction_index,thread_num);
		n2 = mat.refraction_index;
	    }
	}
	Line ref = ray;
	ref.reflect(p,normal);
	ref.direction = -1*ref.direction;
	path->add(p,normal,ref.direction,mat,reflect);
	add_to_path(forward, remaining - change, n2, thread_num, path); 
    }
    else {
	path->add(p,Vector(),Vector(),mat, true);
	return;
    }
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
		Path path(camera+lens_shift);
		add_to_path(ray,ray_depth,1,num,&path);
		c = c+ path.trace();
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
    num_rays = new long[num_threads]();
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
