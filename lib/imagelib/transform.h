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

#ifndef __transform_h__
#define __transform_h__

#include "vector.h"

/* A linear transformation matrix in homogeneous coordinates */
typedef struct {
  double T[3][3];
} trans2D_t;

/* Create a new transform filled with the given values */
trans2D_t *new_transform(double t00, double t01, double t02,
                         double t10, double t11, double t12,
                         double t20, double t21, double t22);

/* Return a new zero transform */
trans2D_t *new_zero_transform();

/* Return a new identity transform */
trans2D_t *new_identity_transform();

/* Create a new 2D translation transform */
trans2D_t *new_translation_transform(double tx, double ty);

/* Create a new 2D rotation transform */
trans2D_t *new_rotation_transform(double theta);

/* Create a new 2D scaling transform */
trans2D_t *new_scaling_transform(double sx, double sy);

/* Create a new 2D transform from the given 3x3 matrix */
trans2D_t *new_transform_vector(double *V);

/* Compute the inverse of the given transform (if it exists) */
trans2D_t *transform_invert(trans2D_t *T);

/* Print out the given transform */
void print_transform(trans2D_t *T);

/* Compute the image of (x,y,1.0) under transformation T, and put the
 * result in (newx, newy) */
void transform_point(trans2D_t *T, double x, double y, double *newx, double *newy);

/* Compute the image of v under T */
v2_t transform_vector(trans2D_t *T, v2_t v);

/* Compute the product of two transforms */
trans2D_t *transform_product(trans2D_t *T, trans2D_t *U);

/* Return a copy of the given transform */
trans2D_t *transform_copy(trans2D_t *T);

/* Free the given transform */
void transform_free(trans2D_t *T);

#endif /* __transform_h__ */
