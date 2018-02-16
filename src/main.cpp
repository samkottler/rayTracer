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

vector<Solid*> objs;
int num_objs = 0;
int num_lights = 0;
default_random_engine* generators;
long* num_rays;
double exposure;
int num_blurs;
int paths_per_pixel;
int width;
int height;
double scale;
double lens_radius;
double focal_length;
int max_ray_depth;
int num_threads;

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

//returns a random ray in the hemisphere around point p with with some normal
//Generates a ray in hemisphere around <0,0,1> and rotates it to proper hemisphere
Line get_reflection_ray(const Line& ray, const Point& p, const Vector& normal, double gloss, int thread_num){
    if(generators[thread_num]() < gloss*generators[thread_num].max()){
	Line ref = ray;
	ref.reflect(p,normal);
	ref.direction = ref.direction*-1;
	return ref;
    }
    else{
	Line new_ray(p,Vector());
	double theta = acos(normal[2]); //angle between normal and <0,0,1>
	Vector n;
	n[0]=0;n[1]=0;n[2]=1;
	n=n.cross(normal); // n is the axis of rotation
	n.normalize();
	double sint = sin(theta);
	double cost = cos(theta);
	Vector v;
	//the first random variable should be between 0-1 to avoid picking a cosine weighted vector
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
	    //rotate the vector to be in correct hemisphere
	    new_ray.direction = cost*(v-n*(n.dot(v))) + n*(n.dot(v)) + sint*n.cross(v);
	}
	return new_ray;
    }
}			      

//returns a perfectly specular refraction ray
Line get_refraction_ray(const Line& ray, const Point& p, const Vector& normal, const Material& mat, double refraction_index, int thread_num){
    double n2;
    Vector actual_normal = normal;
    if (fabs(mat.refraction_index-refraction_index) < 0.0001){ // leaving object
	n2 = 1;
	actual_normal = -1*normal;
    }
    else // entering object
	n2 = mat.refraction_index;
    Line ref = ray;
    double nv = actual_normal.dot(ray.direction);
    //check for total internal reflection
    if((n2*n2)/(refraction_index*refraction_index) < 1-nv*nv) ref.reflect(p,actual_normal);
    else ref.refract(p,actual_normal,refraction_index,n2);
    return ref;
}

//recursive path building function
//at each call calculate the next intersection and reflect or refract until limit is reached or path hits a light
void add_to_path(const Line& ray, int remaining, double refraction_index, int thread_num, Path* path){
    num_rays[thread_num]++;
    Point p(INFINITY,INFINITY,INFINITY);
    Material mat = {Color(),0,0,Color(),Color(),0};
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
        Line forward;
	int change = 1;
	double n2 = 1;
	bool reflect = false;
	if (fabs(mat.refraction_index-refraction_index) < 0.0001){ // leaving
	    forward = get_refraction_ray(ray,p,normal,mat,refraction_index,thread_num);
	    change = 0;
	}
	else{ // entering
	    //if object isn't at all transparent, don't refract
	    if (mat.refraction_index == 0 || generators[thread_num]()<generators[thread_num].max()/2){
		forward = get_reflection_ray(ray,p,normal,mat.gloss,thread_num);
		reflect = true;
	    }
	    else{
		forward = get_refraction_ray(ray,p,normal,mat,refraction_index,thread_num);
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

//adds a bloom effect to the image by aplying several gaussian blurs to any pixels that are over exposed
void blur(Color* colors, int radius, double stddev){
    //construct gaussian kernel
    double kernel[(2*radius+1)*(2*radius+1)];
    double denom = 2*stddev*stddev;
    for (int y = -radius; y<=radius; y++){
	for (int x = -radius; x<=radius; x++){
	    kernel[(y+radius)*(2*radius+1)+x+radius] = exp(-(x*x+y*y)/denom)/denom/M_PI;
	}
    }
    Color* blurs[1+num_blurs];
    blurs[0] = new Color[width*height];
    //start with only the pixels that are over exposed. if a pixel intensity is under the exposure there is no blooming
    for (int i = 0; i<width*height; i++){
	double r,g,b;
	r=g=b=0;
	if (colors[i].r*exposure > 1) r = exposure;
	if (colors[i].g*exposure > 1) g = exposure;
	if (colors[i].b*exposure > 1) b = exposure;
	blurs[0][i] = Color(r,g,b);
    }
    //apply the gaussian blur several times
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
    //add bloom to regular image
    for (int i = 0; i<width*height; i++){
	colors[i] = colors[i] + blurs[num_blurs][i];
    }
    //clean up memory
    for (int i = 0; i<num_blurs+1; i++){
	delete[] blurs[i];
    }
}

//convert the array of Colors to the final array of rgb int values
void expose(int* img, Color* colors){
    blur(colors, 20, 3);
    for (int y=0; y<height; y++){
	for (int x=0; x<width; x++){
	    img[y*width+x]=colors[y*width+x].to_int(1.0/exposure);
	}
    }
}

//construct a path and then trace back the intensities
int completed = 0; //keep track of completed pixels to display progres
mutex completed_lock;
void do_rays_i(Color* img, int num){
    int local_count = 0; //in case this thread can't aquire the lock, keep going and try again later
    //each thread does an entire row
    for (int y = num; y < height; y+=num_threads){
	for (int x = 0; x < width; x++) {
	    Color c;
	    for (int i = 0; i<paths_per_pixel; i++){
		//get a random point inside the pixel
		double xShift = (double)generators[num]()/generators[num].max()-0.5;
		double yShift = (double)generators[num]()/generators[num].max()-0.5;
		Point pixel((double)(x-width/2+xShift)/scale,(double)(height/2-y+yShift)/scale,0);
		//get a random vector inside the lens
		double theta = (double)generators[num]()/generators[num].max()*2*M_PI;
		double rad = (double)generators[num]()/generators[num].max()*lens_radius;
		Vector lens_shift;
		lens_shift[0] = rad*cos(theta); lens_shift[1] = rad*sin(theta); lens_shift[2] = 0;
		Vector D = pixel-camera; //vector to focal point
		D.normalize();
		Point C = camera+focal_length*D; //focal point
		//ray from lens point to focal point
		Line ray = Line(camera+lens_shift, C);
		//build and trace path
		Path path(camera+lens_shift);
		add_to_path(ray,max_ray_depth,1,num,&path);
		c = c+ path.trace();
	    }
	    //average colors and add to array
	    c = c/paths_per_pixel;
	    img[y*width+x]=c;
	}
	//try to get lock to update number of pixels completed
	if (completed_lock.try_lock()){
	    completed += local_count +1;
	    completed_lock.unlock();
	    local_count = 0;
	}
	else{
	    local_count++;
	}
	//no need for a lock here because it doesn't matter if another thread writes to completed while printing
	if (num==0) cout << "\r" << "Progress: " << fixed << setprecision(2) << ((double)(completed)/height)*100 << "%" << flush;
    }
}

int main(int argc, char** argv){
    objs=*read_json_scene("scene.json");
    cout << "Scale:           " << scale << endl;
    cout << "Width:           " << width << endl;
    cout << "Height:          " << height << endl;
    cout << "Paths per pixel: " << paths_per_pixel << endl;
    cout << "Threads:         " << num_threads << endl;
    cout << "Max ray depth:   " << max_ray_depth << endl;
    cout << "Max rays:        " << (double)width*height*paths_per_pixel*max_ray_depth << endl;
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
    delete[] num_rays;
    delete[] generators;
    delete[] img;
    delete[] colors;
    return 0;
}
