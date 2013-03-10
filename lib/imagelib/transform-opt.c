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

/* transform-opt.c */
/* Methods for finding optimal transforms */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "affine.h"
#include "homography.h"
#include "horn.h"
#include "matrix.h"
#include "transform-opt.h"
#include "vector.h"

/* Find the optimal transform of a given class, given point
 * correspondences */
void find_optimal_transform(transform_class_t tclass, int n,
			    v3_t *r_pts, v3_t *l_pts, double *Tout)
{
    int i;

    switch (tclass) {
	case TRANSFORM_TRANSLATE: {
	    v3_t r_centroid = v3_new(0.0, 0.0, 0.0);
	    v3_t l_centroid = v3_new(0.0, 0.0, 0.0);
	    v3_t diff;

	    /* Just subtract the centroid of r_pts from the centroid of
	     * l_pts */

	    for (i = 0; i < n; i++) {
		r_centroid = v3_add(r_centroid, r_pts[i]);
		l_centroid = v3_add(l_centroid, l_pts[i]);
	    }
	    
	    diff = v3_sub(r_centroid, l_centroid);

	    printf("diff: %0.3f, %0.3f\n", Vx(diff) / n, Vy(diff) / n);

	    Tout[0] = 1.0;  Tout[1] = 0.0;  Tout[2] = Vx(diff) / n;
	    Tout[3] = 0.0;  Tout[4] = 1.0;  Tout[5] = Vy(diff) / n;
	    Tout[6] = 0.0;  Tout[7] = 0.0;  Tout[8] = 1.0;

	    break;
	}
	    
	case TRANSFORM_TRANSLATE_ROTATE: {
	    double *weight = (double *) malloc(sizeof(double) * n);
	    double R[9], T[9];
	    double scale;
	    
	    for (i = 0; i < n; i++)
		weight[i] = 1.0;

	    align_horn(n, r_pts, l_pts, R, T, Tout, &scale, weight);

	    free(weight);
	    break;
	}

	case TRANSFORM_AFFINE:
	    align_affine(n, r_pts, l_pts, Tout);
	    break;
	    
	case TRANSFORM_HOMOGRAPHY:
	    align_homography(n, r_pts, l_pts, Tout, 0);
	    break;

	default:
	    printf("[find_optimal_transform] "
		   "Error: invalid transform class!\n");
	    break;
    }
}
