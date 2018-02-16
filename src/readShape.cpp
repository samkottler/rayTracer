#include <vector>
#include <string>
#include <fstream>
#include <cmath>

using namespace std;

#include "json.hpp"
#include "geometry.hpp"

using json = nlohmann::json;

extern int num_lights;
extern int num_objs;
extern double exposure;
extern int num_blurs;
extern int paths_per_pixel;
extern int width;
extern int height;
extern double scale;
extern double lens_radius;
extern double focal_length;
extern int max_ray_depth;
extern int num_threads;

vector<Solid*>* read_json_scene(string filename){
    vector<Solid*>* objs = new vector<Solid*>;
    ifstream scene(filename);
    json details;
    scene>>details;
    exposure = details["exposure"];
    num_blurs = details["num_blurs"];
    num_lights = 0;
    num_objs = details["num_objs"];
    paths_per_pixel = details["paths_per_pixel"];
    width = details["width"];
    height = details["height"];
    scale = details["scale"];
    lens_radius = details["lens_radius"];
    focal_length = details["focal_length"];
    max_ray_depth = details["max_ray_depth"];
    num_threads = details["num_threads"];
    for(int i = 0; i< num_objs; i++){
	json obj;
	scene>>obj;
	if (obj["is_light"].get<bool>()) num_lights+=1;
	if (obj["shape"].get<string>().compare("sphere")==0){
	    json p = obj["center"];
	    Sphere* s = new Sphere(Point(p[0],p[1],p[2]),obj["radius"]);
	    if (obj["is_light"]){
		json c = obj["color"];
		s->material = {Color(c[0],c[1],c[2]), true, 0,Color(), Color(), 0, Color(), Color(), 0,0};
		s->color = Color(c[0],c[1],c[2]);
		s->is_light = true;
	    }
	    else{
		s->is_light = false;
		json diffuse = obj["diffuse"];
		json specular = obj["specular"];
		json rdiffuse = obj["refraction_diffuse"];
		json rspecular = obj["refraction_specular"];
		s->material = {Color(),false,obj["gloss"], Color(diffuse[0], diffuse[1], diffuse[2]), Color(specular[0], specular[1], specular[2]), obj["specular_exp"], Color(rdiffuse[0], rdiffuse[1], rdiffuse[2]), Color(rspecular[0], rspecular[1], rspecular[2]), obj["refraction_specular_exp"], obj["refraction_index"]};
	    }
	    objs->push_back(s);
	}
	else if (obj["shape"].get<string>().compare("plane")==0){
	    json p = obj["point"];
	    json n = obj["normal"];
	    Vector norm;
	    norm[0] = n[0]; norm[1] = n[1]; norm[2] = n[2];
	    norm.normalize();
	    Plane* pl = new Plane(Point(p[0],p[1],p[2]),norm);
	    if (obj["is_light"]) pl->material1 = pl->material2 = {Color(),true,0,Color(),Color(),0, Color(),Color(), 0,0};
	    else{
		json diffuse = obj["diffuse1"];
		json specular = obj["specular1"];
		pl->is_light = false;
		pl->material1 = {Color(),false, obj["gloss"],Color(diffuse[0], diffuse[1], diffuse[2]), Color(specular[0], specular[1], specular[2]), obj["specular_exp"], Color(), Color(), 0, 0};
		diffuse = obj["diffuse2"];
		specular = obj["specular2"];
		pl->material2 = {Color(),false, obj["gloss"],Color(diffuse[0], diffuse[1], diffuse[2]), Color(specular[0], specular[1], specular[2]), obj["specular_exp"], Color(), Color(), 0, 0};
	    }
	    objs->push_back(pl);
	}
	else{
	    json n = obj["normal"];
	    Vector norm;
	    norm[0] = n[0]; norm[1] = n[1]; norm[2] = n[2];
	    norm.normalize();
	    json verts = obj["verts"];
	    int num = verts.size();
	    Point* p_verts = new Point[num];
	    for (int i = 0; i< num; i++){
		json point = verts[i];
		p_verts[i] = Point(point[0],point[1],point[2]);
	    }
	    Face* f = new Face(norm,num,p_verts);
	    if (obj["is_light"]){
		json c = obj["color"];
		f->is_light = true;
		f->color = Color(c[0],c[1],c[2]);
		f->material = {Color(c[0],c[1],c[2]), true,0, Color(), Color(), 0, Color(), Color(), 0,0};
	    }
	    else{
		json diffuse = obj["diffuse"];
		json specular = obj["specular"];
		f->is_light = false;
		f->material = {Color(),false, obj["gloss"],Color(diffuse[0], diffuse[1], diffuse[2]), Color(specular[0], specular[1], specular[2]), obj["specular_exp"], Color(), Color(), 0, 0};
	    }
	    objs->push_back(f);
	}
    }
    return objs;
}
