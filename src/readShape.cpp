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
extern int pixel_samples;
extern int width;
extern int height;
extern int scale;

vector<Solid*>* read_json_scene(string filename){
    vector<Solid*>* objs = new vector<Solid*>;
    ifstream scene(filename);
    json details;
    scene>>details;
    json a = details["ambient"];
    exposure = details["exposure"];
    double r = a[0];
    double g = a[1];
    double b = a[2];
    r = (1.0/(1-r)-1)/exposure;
    g = (1.0/(1-g)-1)/exposure;
    b = (1.0/(1-b)-1)/exposure;
    ambient = Color(r,g,b);
    num_lights = 0;
    num_objs = details["num_objs"];
    pixel_samples = details["pixel_samples"];
    //exposure = details["exposure"];
    width = details["width"];
    height = details["height"];
    scale = details["scale"];
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
	    Vector<3> norm;
	    norm[0] = n[0]; norm[1] = n[1]; norm[2] = n[2];
	    norm.normalize();
	    Plane* pl = new Plane(Point(p[0],p[1],p[2]),norm);
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
