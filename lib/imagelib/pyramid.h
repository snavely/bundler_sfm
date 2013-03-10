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

/* pyramid.h */
/* Routines for creating a gaussian pyramid from an image */

#ifndef __pyramid_h__
#define __pyramid_h__

#ifdef __cplusplus
extern "C" {
#endif

#include <sys/types.h>

#include "image.h"
#include "vector.h"

typedef struct {
    u_int16_t num_levels;  /* Number of levels of the pyramid */
    u_int16_t w, h;        /* Width and height of the bottom of the pyramid */
    img_t **imgs;          /* Images in the pyramid */
} img_pyr_t;

img_pyr_t *img_pyramid_new(int w, int h);

/* Create a gaussian pyramid from an image, with a number of extra levels */
img_pyr_t *img_create_gaussian_pyramid(img_t *img, int extra);

/* Use linear interpolation to compute the color of the point (x,y) 
 * in the given image pyramid */
fcolor_t img_pyramid_pixel_lerp(img_pyr_t *pyr, double x, double y, 
                                double scale);

#ifdef __cplusplus
}
#endif

#endif /* __pyramid_h__ */
