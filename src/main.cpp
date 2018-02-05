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
int num_blurs;
int pixel_samples;
int width,height,scale;
double lens_radius;
double focal_length;

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
	    Line ref = ray;
	    ref.reflect(p,normal);
	    ref.direction=ref.direction*-1;
	    double limit = 0.8*(double)generators[thread_num]()/generators[thread_num].max();
	    for(int i = 0; i< num; i++){
		Vector<3> v;
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
	c = c*(ref_color+get_direct_diffuse(ray,p,normal,thread_num));// + c*(ambient/8);
    }
    return {c,p};
}

void blur(Color* colors, int radius, double stddev){
    double kernel[(2*radius+1)*(2*radius+1)];
    double denom = 2*stddev*stddev;
    for (int y = -radius; y<=radius; y++){
	for (int x = -radius; x<=radius; x++){
	    kernel[(y+radius)*(2*radius+1)+x+radius] = exp(-(x*x+y*y)/denom)/denom/M_PI;
	    //cout << kernel[(y+radius)*(2*radius+1)+x+radius] << " ";
	}
	//cout << endl;
    }
    //int num_blurs = 2;
    Color* blurs[1+num_blurs];
    blurs[0] = new Color[width*height];
    for (int i = 0; i<width*height; i++){
	double r,g,b;
	r=g=b=0;
	if (colors[i].r*exposure > 1) r = exposure;//colors[i].r - 1.0/exposure;
	if (colors[i].g*exposure > 1) g = exposure;//colors[i].g - 1.0/exposure;
	if (colors[i].b*exposure > 1) b = exposure;//colors[i].b - 1.0/exposure;
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
		//cout<< c.b << endl;
		blurs[i][y*width+x] = c;
	    }
	}
    }
    for (int i = 0; i<width*height; i++){
	colors[i] = colors[i] + blurs[num_blurs][i];
    }
    for (int i = 0; i< 1+num_blurs; i++){
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
    for (int y = num; y < height; y+=NUM_THREADS){
	for (int x = 0; x < width; x++) {
	    //int r=0, g=0, b=0;
	    Color c;
	    for (int i = 0; i<pixel_samples; i++){
		double xShift = (double)generators[num]()/generators[num].max()-0.5;
		double yShift = (double)generators[num]()/generators[num].max()-0.5;
		Point pixel((double)(x-width/2+xShift)/scale,(double)(height/2-y+yShift)/scale,0);
		double theta = (double)generators[num]()/generators[num].max()*2*M_PI;
		double rad = (double)generators[num]()/generators[num].max()*lens_radius;
		Vector<3> lens_shift;
		lens_shift[0] = rad*cos(theta); lens_shift[1] = rad*sin(theta); lens_shift[2] = 0;
		Vector<3> D = pixel-camera;
		D.normalize();
		Point C = camera+focal_length*D;
		Line ray = Line(camera+lens_shift, C);
		c = c + trace(ray, DEPTH, num).color;
	    }
	    c = c/pixel_samples;
	    //cout<<c.r<<" "<<c.g<<" " <<c.b<<" " <<hex<< c.to_int(1.0/exposure)<<endl;
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
    cout << "Shadow samples: " << SHADOW_SAMPLES << endl;
    cout << "Pixel samples:  " << pixel_samples << endl;
    cout << "Scatter rays:   " << SCATTER_SAMPLES << endl;
    cout << "Threads:        " << NUM_THREADS << endl;
    cout << "Depth:          " << DEPTH << endl;
    cout << "Max rays:       " << (double)width*height*pixel_samples*DEPTH*pow(SCATTER_SAMPLES,DEPTH) << endl;
    //objs=*read_json_scene("scene.json");
    int* img = new int[width*height];
    Color* colors = new Color[width*height];
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
    cout << "Rays: " <<defaultfloat<< (double)rays << endl;
    expose(img, colors);
    writeImage((char*)"test.png", width, height, img);
    return 0;
}
