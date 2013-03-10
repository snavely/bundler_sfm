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

/* morphology.h */
/* Morphological operators */

#ifndef __morphology_h__
#define __morphology_h__

#include "image.h"

/* Subtract one image from another */
img_t *img_subtract(img_t *img1, img_t *img2);

/* Apply the hit-or-miss operator with the given structural element
 * and 90 degree rotation */
img_t *img_hit_or_miss(img_t *img, int esize, int *elem, int rot90);

/* Dilate an image a given number of times */
img_t *img_dilate(img_t *img, int num_times);

/* Apply thinning to an image */
img_t *img_thin(img_t *img);

#endif /* __morphology_h__ */
