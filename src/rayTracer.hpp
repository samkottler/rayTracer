#ifndef __MY_RAY_TRACER_HPP
#define __MY_RAY_TRACER_HPP

#define SCALE 30
#define WIDTH 400
#define HEIGHT 400
#define SHADOW_SAMPLES 10
#define PIXEL_SAMPLES 50
#define SCATTER_SAMPLES 4
#define DEPTH 10
#define NUM_THREADS 4

extern "C" {
    int writeImage(char* filename, int width, int height, int* image);
}


#endif
