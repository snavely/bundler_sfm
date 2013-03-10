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

#ifndef __transform3D_h__
#define __transform3D_h__

#include "binimage.h"
#include "binvolume.h"

/* A linear transformation matrix in homogeneous coordinates */
typedef struct {
  double T[4][4];
} trans3D_t;

/* Return a new zero transform */
trans3D_t *new_zero_transform3D();

/* Return a new identity transform */
trans3D_t *new_identity_transform3D();

/* Produce a transformed binary volume, given either a 3D or 2D transform */
bvol_t *transform_bin_volume(bvol_t *b, trans3D_t *T);
bvol_t *transform_bin_volume_2D(bvol_t *b, trans2D_t *T);

/* Print out the given transform */
void print_transform3D(trans3D_t *T);

/* Compute the image of (x,y,1.0) under transformation T, and put the
 * result in (newx, newy) */
void transform_point3D(trans3D_t *T, double in[3], double out[3]);

/* Compute the product of two transformations */
trans3D_t *transform_product3D(trans3D_t *T, trans3D_t *U);

/* Reduce T from a 3D transform to a 2D transform.  Assumes that the
 * left-most column of T is the vector [ 0.0, 0.0, 0.0, 1.0 ]^T */
trans2D_t *transform_lower_dimension(trans3D_t *T);

/* Free the given transformation */
void free_transform3D(trans3D_t *T);

#endif /* __transform3D_h__ */
