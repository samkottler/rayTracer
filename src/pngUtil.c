#include <stdlib.h>
#include <png.h>

int writeImage(char* filename, int width, int height, int* image){
    FILE* f = fopen(filename, "wb");
    png_structp wStruct = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_infop iStruct = png_create_info_struct(wStruct);


    png_init_io(wStruct, f);
    png_set_IHDR(wStruct, iStruct, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
		 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    png_write_info(wStruct, iStruct);

    png_bytep row = (png_bytep) malloc(3*width*sizeof(png_byte));
    for (int y = 0; y < height; y++){
	for (int x = 0; x < width; x++){
	    row[x*3+2] = image[y*width+x]&0xff;
	    row[x*3+1] = (image[y*width+x]>>8)&0xff;
	    row[x*3] = (image[y*width+x]>>16)&0xff;
	}
	png_write_row(wStruct, row);
    }
    png_write_end(wStruct, NULL);
}
