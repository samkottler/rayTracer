#include <vector>
#include <string>
#include <fstream>
#include <cmath>

using namespace std;

#include "json.hpp"
#include "geometry.hpp"

using json = nlohmann::json;

extern Color ambient;
extern int num_lights;
extern int num_objs;
extern double exposure;
extern int num_blurs;
extern int pixel_samples;
extern int width;
extern int height;
extern int scale;
extern double lens_radius;
extern double focal_length;
extern int shadow_samples;
extern int ray_depth;
extern int scatter_samples;
extern int num_threads;
extern double distance_scale;

vector<Solid*>* read_json_scene(string filename){
    vector<Solid*>* objs = new vector<Solid*>;
    ifstream scene(filename);
    json details;
    scene>>details;
    json a = details["ambient"];
    exposure = details["exposure"];
    num_blurs = details["num_blurs"];
    double r = a[0];
    double g = a[1];
    double b = a[2];
    //r = (1.0/(1-r)-1)/exposure;
    //g = (1.0/(1-g)-1)/exposure;
    //b = (1.0/(1-b)-1)/exposure;
    ambient = Color(r,g,b);
    num_lights = 0;
    num_objs = details["num_objs"];
    pixel_samples = details["pixel_samples"];
    //exposure = details["exposure"];
    width = details["width"];
    height = details["height"];
    scale = details["scale"];
    lens_radius = details["lens_radius"];
    focal_length = details["focal_length"];
    shadow_samples = details["shadow_samples"];
    ray_depth = details["ray_depth"];
    scatter_samples = details["scatter_samples"];
    num_threads = details["num_threads"];
    distance_scale = details["distance_scale"];
    for(int i = 0; i< num_objs; i++){
	json obj;
	scene>>obj;
	if (obj["is_light"].get<bool>()) num_lights+=1;
	if (obj["shape"].get<string>().compare("sphere")==0){
	    json p = obj["center"];
	    Sphere* s = new Sphere(Point(p[0],p[1],p[2]),obj["radius"]);
	    json c = obj["color"];
	    s->color = Color(c[0],c[1],c[2]);
	    if (obj["is_light"]) s->material = {Color(), true, 0, 0, 0, 0};
	    else{
		json rc = obj["reflect_color"];
		s->material = {Color(rc[0],rc[1],rc[2]),false, obj["reflection_angle"], obj["diffuse"], obj["specular"], obj["specular_exp"]};
	    }
	    objs->push_back(s);
	}
	else{
	    json p = obj["point"];
	    json n = obj["normal"];
	    Vector norm;
	    norm[0] = n[0]; norm[1] = n[1]; norm[2] = n[2];
	    norm.normalize();
	    Plane* pl = new Plane(Point(p[0],p[1],p[2]),norm);
	    json col1 = obj["color1"];
	    json col2 = obj["color2"];
	    pl->color1 = Color(col1[0],col1[1],col1[2]);
	    pl->color2 = Color(col2[0],col2[1],col2[2]);
	    if (obj["is_light"]) pl->material = {Color(),true,0};
	    else{
		json rc = obj["reflect_color"];
		pl->material = {Color(rc[0],rc[1],rc[2]),false, obj["reflection_angle"], obj["diffuse"], obj["specular"], obj["specular_exp"]};
	    }
	    objs->push_back(pl);
	}
    }
    return objs;
}
