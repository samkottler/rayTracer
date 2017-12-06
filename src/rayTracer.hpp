#ifndef __MY_RAY_TRACER_HPP
#define __MY_RAY_TRACER_HPP

#define SCALE 30
#define WIDTH 400
#define HEIGHT 400
#define SHADOW_SAMPLES 10
#define PIXEL_SAMPLES 50
#define DEPTH 3

extern "C" {
    int writeImage(char* filename, int width, int height, int* image);
}


#endif
