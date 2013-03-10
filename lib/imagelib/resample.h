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

#ifndef __resample_h__
#define __resample_h__

#ifdef __cplusplus
extern "C" {
#endif

#include "bmp.h"
#include "color.h"
#include "image.h"
#include "transform.h"

/* Set the background color for unmapped pixels in resampled images */
void set_resample_background(int r, int g, int b);

/* Use linear interpolation to compute the color of the point (x,y) 
 * in the given image */
fcolor_t pixel_lerp(img_t *img, double x, double y);

/* Use linear interpolation to compute the intensity of the point (x,y) 
 * in the given image */
double pixel_lerp_intensity(img_t *img, double x, double y);

/* Compute the pixel that would be at location (x, y) if the given
 * image were transformed with the inverse of the given transform */
fcolor_t pixel_transform(img_t *img, trans2D_t *Tinv, int x, int y);

/* Create a new image by applying transformation T to img and 
 * resampling */
img_t *img_resample(img_t *img, trans2D_t *T);

/* Correct the radial distortion in an image */
img_t *img_fix_radial_distortion(img_t *img, double k1, double k2, double f);

/* Create a new image by applying transformation T to img and 
 * resampling.  Resize the image so that the whole thing fits when
 * transformed */
img_t *img_resample_bbox(img_t *img, trans2D_t *T);

/* Returns true if the image becomes disconnected under the given
 * homography */
int img_disconnected_under_homography(img_t *img, double *H);

void img_sample_random_pt(img_t *img, trans2D_t *Tinv, img_t *img_old, 
			  int nhood_radius, 
			  int *xout, int *yout, double grad_threshold);

#ifdef __cplusplus
}
#endif

#endif /* __resample_h__ */
