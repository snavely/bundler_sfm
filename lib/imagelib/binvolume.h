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

/* binvolume.h */

#ifndef __binvolume_h__
#define __binvolume_h__

#include <sys/types.h>

#include "bmp.h"
#include "binimage.h"
#include "image.h"

typedef enum {
    BVOL_BITS,
    BVOL_HEIGHTS,
} bvol_type_t;

/* Binary volume */
typedef struct {
    u_int32_t w, h, d;
    float axis_weights[3];
    bvol_type_t type;

    union {
        u_int8_t *bits;
        u_int8_t *heights;
    } data;
} bvol_t;

/* Turn on a bit in the volume */
void bvol_setbit(bvol_t *b, int x, int y, int z);

/* Get the value of a bit in the volume */
u_int8_t bvol_getbit(bvol_t *b, int x, int y, int z);

/* Return the number of points in the volume */
u_int32_t bvol_num_points(bvol_t *b);

/* Create a copy of a binary volume */
bvol_t *bvol_copy(bvol_t *b);

/* Return a real number that is "proportional" to the
 * gradient of the volume at the given point (this only applies to 
 * volumes that are height fields).  The returned value will
 * be in the range [0.0, 1.0] */
double bvol_gradient(bvol_t *b, int x, int y, int z);

/* Create an empty binary volume */
bvol_t *new_bin_volume(int w, int h, int d, bvol_type_t type);

/* Slice a plane out of the given binary volume, storing the result in 
 * bimg */
bimg_t *slice_bin_volume(bvol_t *bvol, int z);

/* Routines for reading / writing a binary volume from a file */
void write_bin_volume(FILE *f, bvol_t *b);
bvol_t *read_bin_volume(FILE *f);

/* Convert a bitmap to a binary volume or vice-versa */
bvol_t *bmp2bvm(bmp_t *bmp);
bmp_t *bvm2bmp(bvol_t *b);

/* Convert an image to a binary volume or vice-versa */
bvol_t *img2bvm(img_t *img);
img_t *bvm2img(bvol_t *b);

/* Free the given binary volume */
void free_bin_volume(bvol_t *b);

#define BVOL_HEIGHT(b,x,y) ((b)->data.heights[(y) * (b)->w + (x)])

#endif /* __binvolume_h__ */
