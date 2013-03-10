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
#include <float.h>

#include "binimage.h"
#include "matrix.h"
#include "transform.h"

static double ident[3][3] = { { 1.0, 0.0, 0.0 },
                              { 0.0, 1.0, 0.0 },
                              { 0.0, 0.0, 1.0 } };

/* Return a new zero transform */
trans2D_t *new_zero_transform() {
    trans2D_t *T = calloc(sizeof(trans2D_t), 1);
    return T;
}

/* Return a new identity transform */
trans2D_t *new_identity_transform() {
    trans2D_t *T = malloc(sizeof(trans2D_t));
    memcpy(T->T, ident, 9 * sizeof(double));
    return T;
}

/* Create a new transform filled with the given values */
trans2D_t *new_transform(double t00, double t01, double t02,
                         double t10, double t11, double t12,
                         double t20, double t21, double t22)
{
    trans2D_t *T = malloc(sizeof(trans2D_t));

    T->T[0][0] = t00; T->T[0][1] = t01; T->T[0][2] = t02;
    T->T[1][0] = t10; T->T[1][1] = t11; T->T[1][2] = t12;
    T->T[2][0] = t20; T->T[2][1] = t21; T->T[2][2] = t22;

    return T;
}

/* Create a new 2D transform from the given 3x3 matrix */
trans2D_t *new_transform_vector(double *V) {
    return new_transform(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], V[8]);
}

trans2D_t *new_translation_transform(double tx, double ty) {
    return new_transform(1.0, 0.0, tx,
                         0.0, 1.0, ty,
                         0.0, 0.0, 1.0);
}

trans2D_t *new_rotation_transform(double theta) {
    return new_transform(cos(theta), -sin(theta), 0.0,
                         sin(theta),  cos(theta), 0.0,
                         0.0,         0.0,        1.0);
}

trans2D_t *new_scaling_transform(double sx, double sy) {
    return new_transform(sx,  0.0, 0.0,
                         0.0, sy,  0.0,
                         0.0, 0.0, 1.0);
}

void transform_free(trans2D_t *T) {
    free(T);
}

#if 0
static double det(trans2D_t *T) {
    /* 00 01 02 
     * 10 11 12
     * 20 21 22 */

#define A (T->T)
    return 
        A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
        A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
        A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
#undef A
}
#endif

trans2D_t *transform_product(trans2D_t *T, trans2D_t *U) {
    int i, j, k;
    trans2D_t *R = new_zero_transform();
    
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                R->T[i][j] += T->T[i][k] * U->T[k][j];
            }
        }
    }

    return R;
}

trans2D_t *transform_invert(trans2D_t *T) {
    trans2D_t *Tinv = new_zero_transform();
    matrix_invert(3, (double *)T->T, (double *)Tinv->T);
    return Tinv;
}

trans2D_t *transform_copy(trans2D_t *T) {
    trans2D_t *Tcopy = new_zero_transform();
    memcpy(Tcopy->T, T->T, sizeof(double) * 9);
    return Tcopy;
}


#if 0
trans2D_t *invert_transform(trans2D_t *T) {
    double d = det(T);
    double cofactor;
    int i, j;
    int one[3] = { 1, 0, 0 };
    int two[3] = { 2, 2, 1 };

    trans2D_t *Tinv = new_identity_transform(), *I;
    
    /* Check if the transform is invertible */
    if (d == 0.0) {
        printf("Singular transformation detected!\n");
        return NULL;
    }
  
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            int parity = (i+j) % 2;
            cofactor = (parity ? -1.0 : 1.0) * 
                (T->T[one[i]][one[j]] * T->T[two[i]][two[j]] - 
                 T->T[one[i]][two[j]] * T->T[two[i]][one[j]]);

            Tinv->T[j][i] = cofactor / d;
        }
    }

    I = transform_product(T, Tinv);
    printf("Tinv:\n");
    print_transform(Tinv);
    printf("I:\n");
    print_transform(I);
    transform_free(I);

    return Tinv;
}
#endif

/* Compute the image of (x,y,1.0) under transformation T, and put the
 * result in (newx, newy) */
void transform_point(trans2D_t *T, double x, double y, double *newx, double *newy) {
    double in[3] = { x, y, 1.0 };
    double out[3], invw;

    out[0] = T->T[0][0] * in[0] + T->T[0][1] * in[1] + T->T[0][2] * in[2];
    out[1] = T->T[1][0] * in[0] + T->T[1][1] * in[1] + T->T[1][2] * in[2];
    out[2] = T->T[2][0] * in[0] + T->T[2][1] * in[1] + T->T[2][2] * in[2];

    if (out[2] != 1.0) {
	invw = 1.0 / out[2];
	*newx = out[0] * invw;
	*newy = out[1] * invw;
    } else {
	*newx = out[0];
	*newy = out[1];
    }
}

/* Compute the image of v under T */
v2_t transform_vector(trans2D_t *T, v2_t v) {
    v2_t out;
    transform_point(T, Vx(v), Vy(v), &(Vx(out)), &(Vy(out)));
    return out;
}

void print_transform(trans2D_t *T) {
    int i,j;
    for (i = 0; i < 3; i++) {
        printf("{");
        for (j = 0; j < 3; j++)
            printf(" %f ", T->T[i][j]);
        printf("}\n");
    }
}

