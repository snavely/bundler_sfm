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

#ifndef __binimage_h__
#define __binimage_h__

#include <sys/types.h>

#include "bmp.h"
#include "image.h"
#include "transform.h"
#include "vector.h"

/* Binary image */
typedef struct {
  u_int32_t w, h;
  u_int8_t *bits;
} bimg_t;

#define DISTSQ(x1,y1,x2,y2) (((x1) - (x2)) * ((x1) - (x2)) + ((y1) - (y2)) * ((y1) - (y2)))

/* Turn on a bit in the image */
void bimg_setbit(bimg_t *b, int x, int y);
/* Get the value of a bit in the image */
u_int8_t bimg_getbit(bimg_t *b, int x, int y);

/* Create an empty binary image */
bimg_t *new_bin_image(int w, int h);
void bimg_clear(bimg_t *b);

/* Write a binary image to a file */
int write_bin_image(FILE *f, bimg_t *b);

/* Read a binary image from a file */
bimg_t *read_bin_image(FILE *f);

/* Free a binary image */
void free_bin_image(bimg_t *b);

/* Open a bmp file and return it as a binary image.  Assumes the input
 * file is monochrome */
bimg_t *bmp2bmg(bmp_t *bmp);
bimg_t *img2bmg(img_t *bmp);

/* Convert a binary image to a bitmap file */
bmp_t *bmg2bmp(bimg_t *b);
img_t *bmg2img(bimg_t *b);

/* Print an ascii representation of the given binary image */
void print_bin_image(bimg_t *b);

/* Produce a transformed binary image */
bimg_t *transform_bin_image(bimg_t *b, trans2D_t *t);

/* Find a linear transformation that aligns image B to image A */
trans2D_t *align_bin_image(bimg_t *a, bimg_t *b);

/* Return a list of all points that are turned on in the binary image*/
iv2_t *bin_image_point_list(bimg_t *b, int *num_pts);

#endif /* __binimage_h__ */
