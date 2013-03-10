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

/* filter.h */
/* Routines for filtering images */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef __filter_h__
#define __filter_h__

#include "image.h"

/* Convolve the image about the given pixel with the give m by n
 * kernel */
double img_pixel_convolve_gs(img_t *img, int x, int y, int m, int n, 
			     double *coeffs, double scale);

/* Convolve the image about the given pixel with the given m by n
 * kernel */
void img_pixel_convolve_rgb(img_t *img, int x, int y, int m, int n, 
			    double *coeffs, double *out, double scale);

/* Filter the given image with the m by n array of filter coefficients */
img_t *img_filter(img_t *img, int m, int n, double *coeffs, double scale);

/* Filter is a 1-D filter */
img_t *img_filter_xy(img_t *img, double *filter, int size, int wrap);

/* Compute the magnitude of the gradient of the image at the point
 * (x,y) using a sobel filter */
void img_gradient_sobel_xy(img_t *img, int x, int y, double *dx, double *dy);
    
/* Compute the magnitude of the gradient of the image at the point
 * (x,y) using a sobel filter.  Store the direction of the gradient
 * in the provided theta pointer. */
double img_gradient_sobel(img_t *img, int x, int y, double *theta);

/* Return the laplacian of an image at a given point */
double img_laplacian(img_t *img, int x, int y);

/* Apply a Gaussian filter with variance sigma to the image */
img_t *img_smooth(img_t *img, double sigma, int wrap);

/* Apply a Gaussian filter with variance sigma to the image */
fimg_t *img_smooth_float(img_t *img, double sigma, int wrap);

/* Compute a filter kernel for a gaussian filter with the *
 * given standard deviation */
double *compute_gaussian_filter(double sigma, double num_dev, int *size);

/* Scale an image */
img_t *img_scale(img_t *img, int scale);
/* Scale an image quickly */
img_t *img_scale_fast(img_t *img, int scale);
    
#ifdef __cplusplus
}
#endif

#endif /* __filter_h__ */
