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

/* pyramid.c */
/* Routines for creating a gaussian pyramid from an image */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "pyramid.h"
#include "resample.h"
#include "util.h"


/* Initialize a new image pyramid with the given width and height */
img_pyr_t *img_pyramid_new(int w, int h) {
    int n, num_levels;
    img_pyr_t *pyr;

    if (!is_power_of_two(w) || w != h) {
	printf("[img_pyramid_new] Error: invalid width and height\n");
	return NULL;
    }
    
    n = w;
    num_levels = ilog2(n);

    /* Initialize the pyramid */
    pyr = (img_pyr_t *)malloc(sizeof(img_pyr_t));
    pyr->w = n;  pyr->h = n;
    pyr->num_levels = num_levels;
    pyr->imgs = (img_t **)malloc(sizeof(img_t *) * num_levels);

    return pyr;
}

/* Creates a gaussian pyramid from an image */
img_pyr_t *img_create_gaussian_pyramid(img_t *img, int extra) {
    /* Assumes that the image size is NxN where N=2^k */
    int num_levels;
    img_pyr_t *pyr;
    int i, n, x, y;
    color_t *pxls[4];

    img_t *i_new = img_upsize_square_power_of_two(img);

    printf("Creating pyramid for (%d,%d)\n", img->w, img->h);

    pyr = img_pyramid_new(i_new->w, i_new->h);
    num_levels = pyr->num_levels;
    printf("num_levels = %d\n", num_levels);
    n = i_new->w;

    /* Copy over the first image */
    pyr->imgs[extra] = img_copy(i_new);

    /* Resample the bigger images */
    for (i = 0; i < extra; i++) {
	trans2D_t *T = new_scaling_transform((1 << (extra - i)), (1 << (extra - i)));
	pyr->imgs[i] = img_resample_bbox(img, T);
	pyr->imgs[i]->origin = v2_scale((1 << (extra - i)), img->origin);
	transform_free(T);
    }

    /* Resample the smaller images */
    for (i = extra + 1; i < num_levels; i++) {
	pyr->imgs[i] = img_new((n >> (i - extra)), (n >> (i - extra)));
	pyr->imgs[i]->origin = v2_scale(1.0 / (1 << (i - extra)), img->origin);

	for (y = 0; y < (n >> (i - extra)); y++) {
	    for (x = 0; x < (n >> (i - extra)); x++) {
		color_t p;

		int xmin = (x << 1);
		int ymin = (y << 1);

		if (img_region_is_valid(pyr->imgs[i-1], xmin, xmin + 1, 
                                        ymin, ymin + 1)) {
		    img_get_pixel_square(pyr->imgs[i-1], 
                                         x << 1, y << 1, pxls);
		    
		    p.r = (pxls[0]->r + pxls[1]->r + pxls[2]->r + pxls[3]->r) >> 2;
		    p.g = (pxls[0]->g + pxls[1]->g + pxls[2]->g + pxls[3]->g) >> 2;
		    p.b = (pxls[0]->b + pxls[1]->b + pxls[2]->b + pxls[3]->b) >> 2;
		    
		    img_set_pixel(pyr->imgs[i], x, y, p.r, p.g, p.b);
		}
	    }
	}
    }

#if 1
    for (i = 0; i < num_levels; i++) {
	img_t *img_tmp = pyr->imgs[i];
	pyr->imgs[i] = img_shrink_wrap(pyr->imgs[i]);
	img_free(img_tmp);
    }
#endif

    return pyr;
}

/* Use linear interpolation to compute the color of the point (x,y) 
 * in the given image pyramid */
fcolor_t img_pyramid_pixel_lerp(img_pyr_t *pyr, double x, double y, 
                                double scale) 
{
    double log_scale = 0.0;
    u_int32_t scale1, scale2;
    double t, factor;
    fcolor_t col1, col2, col;

    if (scale <= 1.0) {
        return pixel_lerp(pyr->imgs[0], x, y);
    }

    log_scale = log(scale / log(2.0));

    scale1 = iround(floor(log_scale));
    scale2 = scale1 + 1;

    t = log_scale - (double) scale1;

    if (scale2 >= pyr->num_levels) {
        factor = 1.0 / ((2 << pyr->num_levels) - 1);
        return pixel_lerp(pyr->imgs[pyr->num_levels-1], 
                          x * factor, y * factor);
    }

    factor = 1.0 / (2 << scale1);

    col1 = pixel_lerp(pyr->imgs[scale1], x * factor, y * factor);
    col2 = pixel_lerp(pyr->imgs[scale2], 0.5 * x * factor, 0.5 * x * factor);

    col.r = (1.0 - t) * col1.r + t * col1.r;
    col.g = (1.0 - t) * col1.g + t * col1.g;
    col.b = (1.0 - t) * col1.b + t * col1.b;

    return col;
}

