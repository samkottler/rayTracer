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
    RGB to_RGB() const{
	double r =  2.3706743*x - 0.9000405*y - 0.4706338*z;
	double g = -0.5138850*x + 1.4253036*y + 0.0885814*z;
	double b =  0.0052982*x - 0.0146949*y + 1.0093968*z;
	return RGB(r,g,b);
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
	    ss.ignore(); ss.ignore();
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
};
