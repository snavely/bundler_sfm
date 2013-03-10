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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>

// #include "minpack.h"

#include "binvolume.h"
#include "dist.h"
#include "distutil.h"
#include "lerp.h"
#include "matrix.h"
#include "metric.h"
#include "transform3D.h"
#include "triangle.h"
#include "util.h"

static double ident[4][4] = { { 1.0, 0.0, 0.0, 0.0 },
                              { 0.0, 1.0, 0.0, 0.0 },
                              { 0.0, 0.0, 1.0, 0.0 },
                              { 0.0, 0.0, 0.0, 1.0 } };

/* Return a new zero transform */
trans3D_t *new_zero_transform3D() {
    trans3D_t *T = calloc(sizeof(trans3D_t), 1);
    return T;
}

/* Return a new identity transform */
trans3D_t *new_identity_transform3D() {
    trans3D_t *T = malloc(sizeof(trans3D_t));
    memcpy(T->T, ident, 16 * sizeof(double));
    return T;
}

void free_transform3D(trans3D_t *T) {
    free(T);
}

trans3D_t *transform_product3D(trans3D_t *T, trans3D_t *U) {
    int i, j, k;
    trans3D_t *R = new_zero_transform3D();
    
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            for (k = 0; k < 4; k++) {
                R->T[i][j] += T->T[i][k] * U->T[k][j];
            }
        }
    }

    return R;
}

trans3D_t *invert_transform3D(trans3D_t *T) {
    trans3D_t *Tinv = new_zero_transform3D();
    matrix_invert(4, (double *)T->T, (double *)Tinv->T);
    return Tinv;
}

/* Compute the image of (x,y,1.0) under transformation T, and put the
 * result in (newx, newy) */
void transform_point3D(trans3D_t *T, double in[3], double out[3]) {
    double tmp[4];
    double in_homo[4] = { in[0], in[1], in[2], 1.0 };
    
    int i, j;

    for (i = 0; i < 4; i++) {
        tmp[i] = 0.0;
        for (j = 0; j < 4; j++) {
            tmp[i] += T->T[i][j] * in_homo[j];
        }
    }

    out[0] = tmp[0] / tmp[3];
    out[1] = tmp[1] / tmp[3];
    out[2] = tmp[2] / tmp[3];
}

bvol_t *transform_bin_volume(bvol_t *b, trans3D_t *T) {
    /* Test each point in the plane for inclusion in the new
     * (transformed) binary image */
    bvol_t *Tb = new_bin_volume(b->w, b->h, b->d, BVOL_BITS);

    /* Compute the inverse of T */
    trans3D_t *Tinv = invert_transform3D(T);

    int x, y, z;
    int xc, yc, zc;

    for (z = 0; z < (int) b->d; z++) {
        for (y = 0; y < (int) b->h; y++) {
            for (x = 0; x < (int) b->w; x++) {
                double in[3] = { (double) x, (double) y, (double) z };
                double out[3];
                
                transform_point3D(Tinv, in, out);

                /* Check if the point lies in the set */
                xc = iround(out[0]);
                yc = iround(out[1]);
                zc = iround(out[2]);

                if (bvol_getbit(b, xc, yc, zc))
                    bvol_setbit(Tb, x, y, z);
            }
        }
    }
    
    free_transform3D(Tinv);

    return Tb;
}

bvol_t *transform_bin_volume_2D(bvol_t *b, trans2D_t *T) {
    /* Test each point in the plane for inclusion in the new
     * (transformed) binary image */
    bvol_t *Tb = new_bin_volume(b->w, b->h, b->d, BVOL_BITS);

    /* Compute the inverse of T */
    trans2D_t *Tinv = transform_invert(T);

    int x, y, z;
    int xc, yc, zc;

    for (z = 0; z < (int) b->d; z++) {
        for (y = 0; y < (int) b->h; y++) {
            for (x = 0; x < (int) b->w; x++) {
                double out[3];
                
                transform_point(Tinv, (double) x, (double) y, out + 0, out + 1);
		out[2] = z;

                /* Check if the point lies in the set */
                xc = iround(out[0]);
                yc = iround(out[1]);
                zc = iround(out[2]);

                if (bvol_getbit(b, xc, yc, zc))
                    bvol_setbit(Tb, x, y, z);
            }
        }
    }
    
    transform_free(Tinv);

    return Tb;
}

/* Reduce T from a 3D transform to a 2D transform.  Assumes that the
 * left-most column of T is the vector [ 0.0, 0.0, 0.0, 1.0 ]^T */
trans2D_t *transform_lower_dimension(trans3D_t *T) {
    trans2D_t *T2D = new_zero_transform();
    int i, j;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            T2D->T[i][j] = T->T[i][j];
    
    return T2D;
}
 
void print_transform3D(trans3D_t *T) {
    int i,j;
    for (i = 0; i < 4; i++) {
        printf("{");
        for (j = 0; j < 4; j++)
            printf(" %f ", T->T[i][j]);
        printf("}\n");
    }
}
