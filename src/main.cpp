#include <iostream>

//#include "Vector.hpp"
#include "rayTracer.hpp"

using namespace std;

int main(int argc, char** argv){
    int img[256*256];
    for (int y = 0; y < 256; y++){
	for (int x = 0; x < 256; x++) {
	    img[y*256+x] = (y<<8)+x;
	}
    }
    writeImage((char*)"test.png", 256, 256, img);
    return 0;
}
