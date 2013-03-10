/* 
 *  Copyright (c) 2008  Noah Snavely (snavely (at) cs.washington.edu)
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

/* image.h */
/* Image routines */

#ifndef __image_h__
#define __image_h__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdarg.h>
#include <sys/types.h>

#include "bmp.h"
#include "color.h"
#include "transform.h"
#include "vector.h"

typedef struct {
    /* Width, height */
    u_int16_t w, h;
  
    /* Pixel data */
    color_t *pixels;

    /* Mask of valid pixels */
    u_int8_t *pixel_mask;

    /* The location of the image origin in world space */
    v2_t origin;
} img_t;

typedef struct {
    /* Width, height */
    u_int16_t w, h;
  
    /* Pixel data */
    float *pixels;

    /* Mask of valid pixels */
    u_int8_t *pixel_mask;

    /* The location of the image origin in world space */
    v2_t origin;
} fimg_t;

/* Create a new image with the given width and height */
img_t *img_new(int w, int h);

/* Create a new float image with the given width and height */
fimg_t *fimg_new(int w, int h);

/* Get/set a pixel in the image to the given color */
color_t img_get_pixel(img_t *img, int x, int y);
void img_set_pixel(img_t *img, int x, int y, 
                   u_int8_t r, u_int8_t g, u_int8_t b);
void img_set_RGBA(img_t *img, int x, int y, 
		  u_int8_t r, u_int8_t g, u_int8_t b, u_int8_t a);

float fimg_get_pixel(fimg_t *img, int x, int y);
void fimg_set_pixel(fimg_t *img, int x, int y, float p);

/* Set the alpha on a pixel */
void img_set_alpha(img_t *img, int x, int y, u_int8_t a);

/* Fill in the given array with the 2x2 square of pixels whose
 * bottom-left corner is (x,y) */
void img_get_pixel_square(img_t *img, int x, int y, color_t **sq);

int img_pixel_in_range(img_t *img, double x, double y);

fcolor_t img_average_2x2_square(img_t *img, int x, int y);

/* Set a flag on the given pixel indicating that it should be left
 * out of the image */
// void img_nullify_pixel(img_t *img, int x, int y);

void img_set_valid_pixel(img_t *img, int x, int y);
void fimg_set_valid_pixel(fimg_t *img, int x, int y);

/* Returns 1 is the given pixel is not considered part of the image */
int img_pixel_is_valid(img_t *img, int x, int y);

/* Invalidate the given pixel */
int img_invalidate_pixel(img_t *img, int x, int y);

/* Make a deep copy of the given image */
img_t *img_copy(img_t *img);

/* Add white noise with the given amplitude to the image.
 * The amplitude may be between 0.0 and 1.0 */
void img_add_white_noise(img_t *img, double a);

/* Create gradient maps for the given image */
void img_create_gradient_maps(img_t *img, img_t **x_out, img_t **y_out);

/* Create a subimage of the given image beginning at the 
 * given x, y offset and with width, height equal to w, h */
img_t *img_sub_image(img_t *img, int x, int y, int w, int h);

/* Pad an image with empty pixels */
img_t *img_pad_top(img_t *img, int pixels);
img_t *img_pad_bottom(img_t *img, int pixels);
img_t *img_pad_left(img_t *img, int pixels);
img_t *img_pad_right(img_t *img, int pixels);

/* Merge the given images into a single image stretching horizontally */
img_t *img_merge(img_t *i1, img_t *i2);

/* Fuse the list of images into a single image.  
 * num = number of arguments */
img_t *img_fuse(int num, int blend, ...);

/* Paste the source image into the destination image (in place) */
void img_paste(img_t *dest, img_t *src);
/* Add the second image to the first */
void img_blend(img_t *dest, img_t *src);

/* Compute the magnitude of the gradient of the image at the point (x,y) */
double img_gradient(img_t *img, int x, int y);

/* Normalize an image, returning the result as a new image */
img_t *img_normalize(img_t *img);

/* Equalize an image, returning the result as a new image */
img_t *img_equalize(img_t *img);

/* Compute the mean color of the neighborhood of radius `rad'
 * surrounding the pixel xc, yc */
void img_mean_color(img_t *img, int xc, int yc, int rad, double *mean);

/* Compute the mean grayscale color of the neighborhood of radius
 * `rad' surrounding the pixel xc, yc */
double img_mean_grayscale(img_t *img, int xc, int yc, int rad);

/* Compute the variance of a neighborhood of radius `rad' surrounding
 * the pixel xc, yc */
double img_variance_grayscale(img_t *img, int xc, int yc, int rad);
double img_variance(img_t *img, int xc, int yc, int rad);

/* Find the maximum and minimum variances of all neighborhoods of size
 * rad in the image */
void img_find_min_max_variance(img_t *img, int rad, double *min, double *max);

/* Combine a new color with the existing color of pixel (x,y) in the
 * image, weighted to give equal weight to all parts already
 * contributing to the pixel color */
void pixel_combine(img_t *img, int x, int y, double r, double g, double b, int num_parts);

/* Routines for reading/writing an image to/from a file */
img_t *img_read(FILE *f);
void img_write(FILE *f, img_t *img);
img_t *img_read_bmp_file(char *fname);
void img_write_bmp_file(img_t *img, char *fname);
void img_write_bmp_file_gray(img_t *img, char *fname, int blue);

/* Routines for converting between images and bitmaps */
img_t *bmp2img(bmp_t *bmp);
bmp_t *img2bmp(img_t *img);
bmp_t *img2bmp32(img_t *img);

img_t *fimg2img(fimg_t *img);

/* Convert from an image to grayscale bitmap */
bmp_t *img2bmp_gray(img_t *img, int blue);

/* Convert the image to grayscale (returning a new image) */
img_t *img_convert_grayscale(img_t *img);
    
/* Draw a square on the image */
void img_draw_pt(img_t *img, int x, int y, int size, u_int8_t r, u_int8_t g, u_int8_t b);

/* Draw a line from (x0, y0) to (x1, y1) */
void img_draw_line(img_t *img, int x0, int y0, int x1, int y1, u_int8_t r, u_int8_t g, u_int8_t b);

/* Check if the image of the neighborhood under the given transform
 * belongs to the original image.  This is true if all four corners of
 * the neighborhood are inside the original image */
int nhood_image_within_image(int w, int h, v2_t origin, trans2D_t *T, double x, double y, int nhood_radius);

/* Upsize the image so that is a square and width and height are a
 * power of two */
img_t *img_upsize_square_power_of_two(img_t *img);

/* Upsize the image so that width and height are a power of two */
img_t *img_upsize_power_of_two(img_t *img);

/* Returns true is the given subregion of an image consists of valid
 * pixels */
int img_region_is_valid(img_t *img, int xmin, int xmax, int ymin, int ymax);

/* Light up invalid pixels */
void img_light_invalid_pixels(img_t *img);

/* Turn off all pixels of a given color */
void img_set_invalid_pixels(img_t *img, unsigned char r, unsigned char g, unsigned char b);

/* Compute the width, height, and origin of the image under the given
 * transform */
void img_compute_Tstats(img_t *img, trans2D_t *T, int *w_out, int *h_out, v2_t *origin_out);

/* Free an image */
void img_free(img_t *img);
void fimg_free(fimg_t *img);

/* Convert an image to a 16-bit raw format */
void img_write_16bit_raw(img_t *img, char *filename);

/* Shrink-wrap an image */
img_t *img_shrink_wrap(img_t *img);

/* Copy the image data into a buffer */
void img_buffer(img_t *img, unsigned char *buffer);

/* Copy the image data into a buffer, with alpha */
void img_buffer_alpha(img_t *img, unsigned char *buffer, unsigned char alpha);

/* Copy the (vertically flipped) image data into a buffer */
void img_buffer_flip(img_t *img, unsigned char *buffer);

/* Convert the given image coordinates to NDC (normalized device
 * coordinates */
void img_to_NDC(img_t *img, double x, double y, double *x_out, double *y_out);

/* Convert the given image coordinates from NDC (normalized device
 * coordinates */
void img_from_NDC(img_t *img, double x, double y, 
		  double *x_out, double *y_out);

/* Convert an image to a locally normalized form */
img_t *img_normalize_local(img_t *img, double sigma);
    
#define IMG_GET_PIXEL(i,x,y) (i)->pixels[(y) * (i)->w + (x)]

#ifdef __cplusplus
}
#endif

#endif /* __image_h__ */
