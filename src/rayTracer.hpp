#ifndef __MY_RAY_TRACER_HPP
#define __MY_RAY_TRACER_HPP

#define SCALE 15
#define WIDTH 200
#define HEIGHT 200
#define SHADOW_SAMPLES 10
#define PIXEL_SAMPLES 50
#define SCATTER_SAMPLES 50
#define DEPTH 3
#define NUM_THREADS 4

extern "C" {
    int writeImage(char* filename, int width, int height, int* image);
}


#endif
