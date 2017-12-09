class RGB{
public:
    double r,g,b;
    RGB(double r1, double g1, double b1){
	r = r1;
	g = g1;
	b = b1;
    }
    int to_int() const{
	return (((int)(r*255))<<16) + (((int)(g*255))<<8) + ((int)(b*255));
    }
};

class XYZ{
public:
    double x,y,z;
    XYZ(double x1, double y1, double z1){
	x = x1;
	y = y1;
	z = z1;
    }
    XYZ(){
	x=y=z=0;
    }
    RGB to_RGB() const{
	double rgb[3];
	rgb[0] =  2.3706743*x - 0.9000405*y - 0.4706338*z;
	rgb[1] = -0.5138850*x + 1.4253036*y + 0.0885814*z;
	rgb[2] =  0.0052982*x - 0.0146949*y + 1.0093968*z;
	for(int i  = 0; i< 3; i++){
	    if (rgb[i]<=0.0031308) rgb[i]*=12.92;
	    else rgb[i] = 1.055*pow(rgb[i],1/2.4)-0.055;
	    if (rgb[i]<0) rgb[i] = 0;
	    if (rgb[i]>1) rgb[i] = 1;
	}
	return RGB(rgb[0],rgb[1],rgb[2]);
    }
    void normalize(double n){
	x/=n;
	y/=n;
	z/=n;
    }
};

class Spectrum{
public:
    double samples[(int)440/SPECTRUM_RESOLUTION];
    static double x[440];
    static double y[440];
    static double z[440];
    static bool initialized;
    void initialize(){
	ifstream file("data/lin2012xyz2e_1_7sf.csv");
	string line;
	for(int i= 0; i<440; i++){
	    getline(file,line);
	    stringstream ss(line);
	    int w;
	    ss>>w;
	    ss.ignore();
	    ss>>x[i]; ss.ignore();
	    ss>>y[i]; ss.ignore();
	    ss>>z[i]; ss.ignore();
	}
	initialized = true;
    }
    XYZ to_XYZ(){
	if (!initialized)
	    initialize();
	double X=0,Y=0,Z=0, total=0;
	for(int i = 0; i<440/SPECTRUM_RESOLUTION; i++){
	    X+=samples[i]*x[i*SPECTRUM_RESOLUTION];
	    Y+=samples[i]*y[i*SPECTRUM_RESOLUTION];
	    Z+=samples[i]*z[i*SPECTRUM_RESOLUTION];
	    total+=y[i*SPECTRUM_RESOLUTION];
	}
	return XYZ(X/total, Y/total, Z/total);
    }
    void operator*= (const Spectrum& other){
	for (int i = 0; i< 440/SPECTRUM_RESOLUTION; i++){
	    samples[i]*=other.samples[i];
	}
    }
    void operator*= (double d){
	for (int i = 0; i< 440/SPECTRUM_RESOLUTION; i++){
	    samples[i]*=d;
	}
    }
    void operator+= (const Spectrum& other){
	for (int i = 0; i< 440/SPECTRUM_RESOLUTION; i++){
	    samples[i]+=other.samples[i];
	}
    }
    void operator/= (double d){
	for (int i = 0; i< 440/SPECTRUM_RESOLUTION; i++){
	    samples[i]/=d;
	}
    }
};
bool Spectrum::initialized = false;
double Spectrum::x[440];
double Spectrum::y[440];
double Spectrum::z[440];
