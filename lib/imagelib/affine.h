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

/* affine.h */
/* Computes an affine transformation */

#ifndef __affine_h__
#define __affine_h__

#ifdef __cplusplus
extern "C" {
#endif

#include "vector.h"

/* Computes the affine transformation that, when applied to the points
 * in l_pts, minimizes the least-squares error between the result and
 * the corresponding points in r_pts.
 * 
 * n -- number of points
 * r_pts -- matches
 * l_pts -- initial points 
 * Tout -- on return, contains the 3x3 transformation matrix */
void align_affine(int num_pts, v3_t *r_pts, v3_t *l_pts, double *Tout);

void align_affine_3D(int num_pts, v3_t *r_pts, v3_t *l_pts, double *Tout);
    
#ifdef __cplusplus
}
#endif

#endif /* __affine_h__ */
