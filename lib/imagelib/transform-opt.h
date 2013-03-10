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

/* transform-opt.h */
/* Methods for finding optimal transforms */

#ifndef __transform_opt_h__
#define __transform_opt_h__

/* List of known transformation classes */
typedef enum {
    TRANSFORM_TRANSLATE,
    TRANSFORM_TRANSLATE_ROTATE,
    TRANSFORM_AFFINE,
    TRANSFORM_HOMOGRAPHY
} transform_class_t;

/* Find the optimal transform of a given class, given point
 * correspondences */
void find_optimal_transform(transform_class_t tclass, int n,
			    v3_t *r_pts, v3_t *l_pts, double *Tout);

#endif /* __transform_opt_h__ */
