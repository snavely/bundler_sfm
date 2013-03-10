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

/* filter.c */
/* Routines for filtering images */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bmp.h"
#include "color.h"
#include "defines.h"
#include "filter.h"
#include "image.h"
#include "matrix.h"
#include "util.h"

/* Convolve the image about the given pixel with the give m by n
 * kernel */
double img_pixel_convolve_gs(img_t *img, int x, int y, int m, int n, 
			     double *coeffs, double scale)
{
    int xmin = x - (n / 2), xmax = x + (n / 2);
    int ymin = y - (m / 2), ymax = y + (m / 2);
    
    int xi, yi;
    double sum = 0.0;

    for (yi = ymin; yi <= ymax; yi++) {
	for (xi = xmin; xi <= xmax; xi++) {
	    int i = yi - ymin, j = xi - xmin;
	    color_t c = img_get_pixel(img, xi, yi);
	    double col = color_intensity(c);
	    double coeff = coeffs[i * n + j];

	    sum += coeff * col;
	}
    }

    return sum / scale;
}


void img_pixel_convolve_rgb(img_t *img, int x, int y, int m, int n, 
			    double *coeffs, double *out, double scale) {
    int xmin = x - (n / 2), xmax = x + (n / 2);
    int ymin = y - (m / 2), ymax = y + (m / 2);
    
    int xi, yi;
    double sum[3] = { 0.0, 0.0, 0.0 };

    if (xmin < 0 || xmax > img->w - 1 || ymin < 0 || ymax > img->h - 1)
        printf("[img_pixel_convolve_rgb] Pixel out of bounds\n");

    for (yi = ymin; yi <= ymax; yi++) {
	for (xi = xmin; xi <= xmax; xi++) {
	    int i = yi - ymin, j = xi - xmin;
	    color_t c = img_get_pixel(img, xi, yi);
	    double coeff = coeffs[i * n + j];

	    sum[0] += coeff * c.r;
	    sum[1] += coeff * c.g;
	    sum[2] += coeff * c.b;
	}
    }

    out[0] = sum[0] / scale;
    out[1] = sum[1] / scale;
    out[2] = sum[2] / scale;
}

/* Filter the given image with the m by n array of filter coefficients */
img_t *img_filter(img_t *img, int m, int n, double *coeffs, double scale) {
    double c[3];
    int x, y;
    img_t *img_out = img_new(img->w, img->h);

    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
            if (x < n/2 || y < m/2 || x >= img->w - n/2 || y >= img->h - m/2)
                img_set_pixel(img_out, x, y, 0, 0, 0);
            else {
                img_pixel_convolve_rgb(img, x, y, m, n, coeffs, c, scale);
                img_set_pixel(img_out, x, y, 
                              iround(clamp(c[0], 0.0, 256.0)), 
                              iround(clamp(c[1], 0.0, 256.0)), 
                              iround(clamp(c[2], 0.0, 256.0)));
            }
	}
    }

    return img_out;
}


static double sobel_x[9] = { -1, 0, 1,
			     -2, 0, 2,
			     -1, 0, 1 };

static double sobel_y[9] = { -1, -2, -1,
 			      0,  0,  0,
			      1,  2,  1 };


/* Compute the magnitude of the gradient of the image at the point
 * (x,y) using a sobel filter */
void img_gradient_sobel_xy(img_t *img, int x, int y, double *dx, double *dy) {
    *dx = img_pixel_convolve_gs(img, x, y, 3, 3, sobel_x, 1.0);
    *dy = img_pixel_convolve_gs(img, x, y, 3, 3, sobel_y, 1.0);
}

/* Compute the magnitude of the gradient of the image at the point
 * (x,y) using a sobel filter */
double img_gradient_sobel(img_t *img, int x, int y, double *theta) {
    double dx = img_pixel_convolve_gs(img, x, y, 3, 3, sobel_x, 1.0);
    double dy = img_pixel_convolve_gs(img, x, y, 3, 3, sobel_y, 1.0);

    if (theta)
	*theta = atan2(dy, dx);

    return sqrt(dx * dx + dy * dy);
}

/* Laplace filter kernel */
static double laplace_kernel[25] = 
  { 0,  0,   1,  0,  0,
    0,  1,   2,  1,  0,
    1,  2, -16,  2,  1,
    0,  1,   2,  1,  0,
    0,  0,   1,  0,  0 };

/* Return the laplacian of an image at a given point */
double img_laplacian(img_t *img, int x, int y) {
  return img_pixel_convolve_gs(img, x, y, 5, 5, laplace_kernel, 1.0);
}

#ifdef WIN32
static double erf (double x) {
  int sign;
  double t;

  if (x < 0)
    {
      sign = -1;
      x = -x;
    }
  else sign = 1;
   
  /* This form uses equation 7.126 of A&S.  It is accurate to 1.5e-7. */
   
  t = 1 / (1 + 0.3275911 * x);
   
  if (x < 1.0e17) {
    t = 1.0 - (t * (0.254829592 + t * 
		    (-0.284496736 + t *
		     (1.421413741 + t *
		      (-1.453152027 + t *
		       (1.061405429)))))) * exp (-x * x);
  }
  else t = 1.0;

  t = t * sign;
   
  return t;
} 
#endif

/* Compute a filter kernel for a gaussian filter with the *
 * given standard deviation */
double *compute_gaussian_filter(double sigma, double num_dev, int *size) {
    /* Go to three standard deviations */
    int rad = (int)ceil(sigma * num_dev), x;
    int my_size = 1 + 2 * rad;

    double *filter = (double *)malloc(sizeof(double) * my_size);
    double sum = 0.0;

    /* Fill in the filter */
    for (x = -rad; x <= rad; x++) {
	filter[x + rad] = 0.5 * (erf((x + 0.5) / (sqrt(2.0) * sigma)) - erf((x - 0.5) / (sqrt(2.0) * sigma)));
	sum += filter[x + rad];
    }

    /* Normalize the filter */
    for (x = 0; x < my_size; x++) 
	filter[x] /= sum;
    
    *size = my_size;

    return filter;
}


/* Filter is a 1-D filter */
img_t *img_filter_xy(img_t *img, double *filter, int size, int wrap) {
    int x, y, i, idx;
    int rad = size / 2;
    int w = img->w, h = img->h;

    double *tmp_img = (double *)malloc(sizeof(double) * 3 * w * h);
    img_t *new_img = img_new(w, h);

    /* Filter in the x direction */
    for (y = 0; y < h; y++) {
	color_t *row = img->pixels + y * w;
	
	for (x = 0; x < w; x++) {
	    double rsum = 0.0, gsum = 0.0, bsum = 0.0;
	    
	    for (i = -rad; i <= rad; i++) {
		if (wrap)
		    idx = (x + i + w) % w;
		else
		    idx = CLAMP(x + i, 0, w - 1);
		
		rsum += filter[i+rad] * ((double) row[idx].r);
		gsum += filter[i+rad] * ((double) row[idx].g);
		bsum += filter[i+rad] * ((double) row[idx].b);
	    }

	    idx = y * 3 * w + 3 * x;

	    tmp_img[idx+0] = rsum;
	    tmp_img[idx+1] = gsum;
	    tmp_img[idx+2] = bsum;
	}
    }

    /* Filter in the y direction */
    for (x = 0; x < w; x++) {
	for (y = 0; y < h; y++) {
	    double rsum = 0.0, gsum = 0.0, bsum = 0.0;
	    
	    for (i = -rad; i <= rad; i++) {
		if (wrap)
		    idx = ((y + i + h) % h) * 3 * w + x * 3;
		else
		    idx = (CLAMP(y + i, 0, h - 1)) * 3 * w + x * 3;

		rsum += filter[i+rad] * ((float) tmp_img[idx + 0]);
		gsum += filter[i+rad] * ((float) tmp_img[idx + 1]);
		bsum += filter[i+rad] * ((float) tmp_img[idx + 2]);
	    }

	    idx = y * 3 * w + 3 * x;
	    
	    img_set_pixel(new_img, x, y, 
			  CLAMP(iround(rsum), 0, 255),
			  CLAMP(iround(gsum), 0, 255),
			  CLAMP(iround(bsum), 0, 255));
	}
    }

    free(tmp_img);

    return new_img;
}

/* Filter is a 1-D filter */
fimg_t *img_filter_xy_float(img_t *img, double *filter, int size, int wrap) {
    int x, y, i, idx;
    int rad = size / 2;
    int w = img->w, h = img->h;

    double *tmp_img = (double *)malloc(sizeof(double) * w * h);
    fimg_t *new_img = fimg_new(w, h);

    /* Filter in the x direction */
    for (y = 0; y < h; y++) {
	color_t *row = img->pixels + y * w;
	
	for (x = 0; x < w; x++) {
	    double sum = 0.0;
	    
	    for (i = -rad; i <= rad; i++) {
		if (wrap)
		    idx = (x + i + w) % w;
		else
		    idx = CLAMP(x + i, 0, w - 1);
		
		sum += filter[i+rad] * ((double) row[idx].r);
	    }

	    idx = y * w + x;

	    tmp_img[idx] = sum;
	}
    }

    /* Filter in the y direction */
    for (x = 0; x < w; x++) {
	for (y = 0; y < h; y++) {
	    double sum = 0.0;
	    
	    for (i = -rad; i <= rad; i++) {
		if (wrap)
		    idx = ((y + i + h) % h) * w + x;
		else
		    idx = (CLAMP(y + i, 0, h - 1)) * w + x;

		sum += filter[i+rad] * ((float) tmp_img[idx]);
	    }

	    idx = y * w + x;
	    
	    fimg_set_pixel(new_img, x, y, (float) sum);
	}
    }

    free(tmp_img);

    return new_img;
}

/* Apply a Gaussian filter with variance sigma to the image */
img_t *img_smooth(img_t *img, double sigma, int wrap) {
    int size;
    double *filter;
    img_t *new_img;
    
    filter = compute_gaussian_filter(sigma, 2.0, &size);    
    // matrix_print(1, size, filter);
    new_img = img_filter_xy(img, filter, size, wrap);

    free(filter);

    return new_img;
}

/* Apply a Gaussian filter with variance sigma to the image */
fimg_t *img_smooth_float(img_t *img, double sigma, int wrap) {
    int size;
    double *filter;
    fimg_t *new_img;
    
    filter = compute_gaussian_filter(sigma, 2.0, &size);
    // matrix_print(1, size, filter);
    new_img = img_filter_xy_float(img, filter, size, wrap);

    free(filter);

    return new_img;
}

/* Scale an image */
img_t *img_scale(img_t *img, int scale) 
{
    int w_out = img->w / scale + ((img->w % 8) != 0);
    int h_out = img->h / scale + ((img->h % 8) != 0);
    img_t *img_out = img_new(w_out, h_out);
    int x, y;

    /* First filter the image */
    img_t *img_filt = img_smooth(img, 1.0 * ilog2(scale), 0);

    for (y = 0; y < img->h; y += scale) {
	for (x = 0; x < img->w; x += scale) {
	    color_t c = img_get_pixel(img_filt, x, y);
	    img_set_pixel(img_out, x / scale, y / scale, c.r, c.g, c.b);
	}
    }

    img_free(img_filt);
    
    return img_out;
}

/* Scale an image quickly */
img_t *img_scale_fast(img_t *img, int scale) 
{
    int w = img->w;
    int h = img->h;
    int w_out = w / scale + ((w % 8) != 0);
    int h_out = h / scale + ((h % 8) != 0);
    img_t *img_out = img_new(w_out, h_out);
    int x, y;

    /* First filter the image */
    for (y = 0; y < h; y += scale) {
	color_t *row = img->pixels + y * w;

	for (x = 0; x < w; x += scale) {
	    int count = 0;
	    int rsum = 0, gsum = 0, bsum = 0;
	    int xi, yi;

	    color_t *c = row + x;
	    unsigned char r, g, b;

	    for (yi = 0; yi < scale; yi++) {
		for (xi = 0; xi < scale; xi++) {
		    if (x + xi >= w || y + yi >= h)
			goto End;
		    
		    rsum += (c+xi)->r;
		    gsum += (c+xi)->g;
		    bsum += (c+xi)->b;

		    count++;
		}

		c += w;
	    }

	End:

	    if (x / scale >= w_out || y / scale >= h_out)
		continue;

	    r = rsum / count;
	    g = gsum / count;
	    b = bsum / count;
	    img_set_pixel(img_out, x / scale, y / scale, r, g, b);
	}
    }

    return img_out;
}
