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

/* morphology.c */
/* Morphological operators */

#include "image.h"
#include "morphology.h"

#define ELEM_DONT_CARE 2

/* Subtract one image from another */
img_t *img_subtract(img_t *img1, img_t *img2) 
{
    int w = img1->w, h = img1->h;
    int x, y;

    img_t *img_out;

    if (img2->w != w || img2->h != h) {
	printf("[img_subtract] Image dimensions do not match\n");
	return NULL;
    }

    img_out = img_new(w, h);

    for (y = 0; y < h; y++) {
	for (x = 0; x < w; x++) {
	    int on1 = (img_get_pixel(img1, x, y).r == 0xff);
	    int on2 = (img_get_pixel(img2, x, y).r == 0xff);

	    if (on1 && !on2) {
		img_set_pixel(img_out, x, y, 0xff, 0xff, 0xff);
	    } else {
		img_set_pixel(img_out, x, y, 0x0, 0x0, 0x0);
	    }
	}
    }
    
    return img_out;
}

static int rot90m[4][4] = 
    { {  1,  0,  0,  1 },
      {  0,  1, -1,  0 },
      { -1,  0,  0, -1 },
      {  0, -1,  1,  0 } };

/* Dilate an image a given number of times */
img_t *img_dilate(img_t *img, int num_times) 
{
    int w = img->w, h = img->h;
    int x, y, i;

    img_t *img_out = img_new(w, h);
    
    for (i = 0; i < num_times; i++) {
	for (y = 0; y < h; y++) {
	    for (x = 0; x < w; x++) {
		int on = 0;
		int dx, dy;

		for (dy = -1; dy <= 1; dy++) {
		    for (dx = -1; dx <= 1; dx++) {
			if (x+dx < 0 || x+dx >= w || y+dy < 0 || y+dy >= h)
			    continue;
		    
			if (img_get_pixel(img, x+dx, y+dy).r > 0x0) {
			    on = 1;
			    goto FoundOne;
			}
		    }
		}

	    FoundOne:
		if (on == 1) {
		    img_set_pixel(img_out, x, y, 0xff, 0xff, 0xff);
		} else {
		    img_set_pixel(img_out, x, y, 0x0, 0x0, 0x0);
		}
	    }
	}

	if (i > 0)
	    img_free(img);
	
	if (i < num_times - 1)
	    img = img_copy(img_out);
    }

    return img_out;
}

/* Apply the hit-or-miss operator with the given structural element
 * and 90 degree rotation */
img_t *img_hit_or_miss(img_t *img, int esize, int *elem, int rot90) 
{
    int rad = esize / 2;
    int w = img->w, h = img->h;
    int x, y;

    img_t *img_out = img_new(w, h);

    rot90 = rot90 & 0x3;

    for (y = 0; y < h; y++) {
	for (x = 0; x < w; x++) {
	    int dx, dy;
	    int match = 1;

	    for (dy = -rad; dy <= rad; dy++) {
		for (dx = -rad; dx <= rad; dx++) {
		    int p[2] = { dx, dy };
		    int *rot90t = rot90m[rot90];
		    int pr[2] = { rot90t[0] * p[0] + rot90t[1] * p[1],
				  rot90t[2] * p[0] + rot90t[3] * p[1] };

		    int on;
		    int e = elem[(pr[1]+rad) * esize + (pr[0]+rad)];

		    if (x+dx < 0 || x+dx >= w || y+dy < 0 || y+dy >= h) {
			on = 0;
		    } else {
			on = (img_get_pixel(img, x+dx, y+dy).r == 0xff);
		    }

		    if (e == ELEM_DONT_CARE) {
			continue;
		    } else {
			if (on != e) {
			    match = 0;
			}
		    }

		    if (match == 0) break;
		}

		if (match == 0) break;
	    }

	    if (match) {
		img_set_pixel(img_out, x, y, 0xff, 0xff, 0xff);
	    } else {
		img_set_pixel(img_out, x, y, 0x0, 0x0, 0x0);
	    }
	}
    }

    return img_out;
}

/* Apply thinning to an image */
img_t *img_thin(img_t *img) {
    int i;

    int elem1[9] = { 0, 0, 0, 
		     2, 1, 2,
		     1, 1, 1 };

    int elem2[9] = { 2, 0, 0,
		     1, 1, 0,
		     2, 1, 2 };

    img_t *img_curr = img_copy(img);
    
    // printf("thinning...\n");

    for (i = 0; i < 4; i++) {
	img_t *bdry = img_hit_or_miss(img, 3, elem1, i);
	img_t *thin = img_subtract(img_curr, bdry);
	img_t *bdry2 = img_hit_or_miss(thin, 3, elem2, i);
	img_t *thin2 = img_subtract(thin, bdry2);

	img_free(img_curr);
	img_free(bdry);
	img_free(thin);
	img_free(bdry2);
	
	img_curr = thin2;
    }

    return img_curr;
}
