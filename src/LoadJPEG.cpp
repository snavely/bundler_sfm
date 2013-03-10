/* 
 *  Copyright (c) 2008-2010  Noah Snavely (snavely (at) cs.cornell.edu)
 *    and the University of Washington
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 */

/* LoadJPEG.cpp */
/* Code for reading a jpeg file */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <jpeglib.h>

#include "image.h"

void GetJPEGDimensions(const char *filename, int &w, int &h)
{
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;
    
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);
    
    FILE *f;
    
    if ((f = fopen(filename, "rb")) == NULL) {
        printf("[GetJPEGDimensions] Error: can't open file %s for reading\n", 
               filename);
        return;
    }
        
    jpeg_stdio_src(&cinfo, f);
    jpeg_read_header(&cinfo, TRUE);

    w = cinfo.image_width;
    h = cinfo.image_height;

    printf("[GetJPEGDimensions] File %s: ( %d , %d )\n", filename, w, h);
    
    jpeg_destroy_decompress(&cinfo);

    fclose(f);
}

/* Note: information on libjpeg can be found here: 
 * http://www.jpegcameras.com/libjpeg/libjpeg-2.html
 */

img_t *LoadJPEG(const char *filename)
{
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;
    
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);
    
    FILE *f;
    
    if ((f = fopen(filename, "rb")) == NULL) {
        printf("[LoadJPEG] Error: can't open file %s for reading\n", filename);
        return NULL;        
    }
        
    jpeg_stdio_src(&cinfo, f);
    jpeg_read_header(&cinfo, TRUE);
    jpeg_start_decompress(&cinfo);

    int w = cinfo.output_width;
    int h = cinfo.output_height;
    int n = cinfo.output_components;

    assert(n == 1 || n == 3);

    img_t *img = img_new(w, h);
    JSAMPROW row = new JSAMPLE[n * w];

    for (int y = 0; y < h; y++) {
        jpeg_read_scanlines(&cinfo, &row, 1);

        for (int x = 0; x < w; x++) {
            if (n == 3) {
                img_set_pixel(img, x, h - y - 1, 
                              row[3 * x + 0], row[3 * x + 1], row[3 * x + 2]);
            } else if (n == 1) {
                img_set_pixel(img, x, h - y - 1, row[x], row[x], row[x]);
            }
        }
    }
    
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);

    delete [] row;

    fclose(f);

    return img;
}

void WriteJPEG(const img_t *img, const char *filename)
{
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    
    
    FILE *outfile;

    if ((outfile = fopen(filename, "wb")) == NULL) {
        printf("[WriteJPEG] can't open file %s for writing\n", filename);
        return;
    }

    jpeg_stdio_dest(&cinfo, outfile);

    cinfo.image_width = img->w;     /* image width and height, in pixels */
    cinfo.image_height = img->h;
    cinfo.input_components = 3;     /* # of color components per pixel */
    cinfo.in_color_space = JCS_RGB; /* colorspace of input image */
    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, 98, TRUE);

    jpeg_start_compress(&cinfo, TRUE);
    
    JSAMPROW row = new JSAMPLE[3 * img->w];
    for (int y = 0; y < img->h; y++) {
        // JSAMPROW row_pointer[1];        /* pointer to a single row */
        int row_stride;                 /* physical row width in buffer */
        row_stride = img->w * 3;        /* JSAMPLEs per row in image_buffer */

        for (int x = 0; x < img->w; x++) {
            color_t c = img_get_pixel((img_t *) img, x, img->h - y - 1);
            row[3 * x + 0] = c.r;
            row[3 * x + 1] = c.g;
            row[3 * x + 2] = c.b;
        }
        
        jpeg_write_scanlines(&cinfo, &row, 1);
    }

    jpeg_finish_compress(&cinfo);
    jpeg_destroy_compress(&cinfo);

    delete [] row;

    fclose(outfile);
}



#if 0
int main()
{
    img_t *img = ReadJPEG("test.jpg");

    img_write_bmp_file(img, "test.bmp");

    return 0;
}
#endif
