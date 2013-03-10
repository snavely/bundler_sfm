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

/* affine.c */
/* Computes an affine transformation */

#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include "vector.h"

/* Computes the affine transformation that, when applied to the points
 * in l_pts, minimizes the least-squares error between the result and
 * the corresponding points in r_pts.
 * 
 * n -- number of points
 * r_pts -- matches
 * l_pts -- initial points 
 * Tout -- on return, contains the 3x3 transformation matrix */
void align_affine(int num_pts, v3_t *r_pts, v3_t *l_pts, double *Tout) {
    int m = num_pts * 2; /* Rows of A */
    int n = 6;           /* Columns of A */
    int nrhs = 1;        /* Columns of X */

    double *A = malloc(sizeof(double) * m * n);    /* Left-hand matrix */
    double *B = malloc(sizeof(double) * m * nrhs); /* Right-hand matrix */
    // double *X = malloc(sizeof(double) * n * nrhs);
    double xVec[6];

    int i, base;
    
    for (i = 0; i < num_pts; i++) {
	base = 2 * i * n;

	/* Row 1 */
	A[base + 0] = (double) Vx(l_pts[i]);
	A[base + 1] = (double) Vy(l_pts[i]);
	A[base + 2] = 1.0;
	A[base + 3] = A[base + 4] = A[base + 5] = 0.0;

	base = (2 * i + 1) * n;
	/* Row 2 */
	A[base + 0] = A[base + 1] = A[base + 2] = 0.0;
	A[base + 3] = (double) Vx(l_pts[i]);
	A[base + 4] = (double) Vy(l_pts[i]);
	A[base + 5] = 1.0;
	
	B[2 * i + 0] = (double) Vx(r_pts[i]);
	B[2 * i + 1] = (double) Vy(r_pts[i]);
    }

    /* Run the driver to dgelsy */
    dgelsy_driver(A, B, xVec, m, n, nrhs);

    Tout[0] = xVec[0];  Tout[1] = xVec[1];  Tout[2] = xVec[2];
    Tout[3] = xVec[3];  Tout[4] = xVec[4];  Tout[5] = xVec[5];
    Tout[6] = 0.0;  Tout[7] = 0.0;  Tout[8] = 1.0;

    free(A);
    free(B);
    // free(X);
}


/* Computes the 3D affine transformation that, when applied to the points
 * in l_pts, minimizes the least-squares error between the result and
 * the corresponding points in r_pts.
 * 
 * n -- number of points
 * r_pts -- matches
 * l_pts -- initial points 
 * Tout -- on return, contains the 4x4 transformation matrix */

void align_affine_3D(int num_pts, v3_t *r_pts, v3_t *l_pts, double *Tout) {
    if (num_pts < 4) {
	printf("[align_affine_3D] System is underconstrained!\n");
	return;
    } else {
	int m = num_pts * 3; /* Rows of A */
	int n = 12;           /* Columns of A */
	int nrhs = 1;        /* Columns of X */

	double *A = malloc(sizeof(double) * m * n);    /* Left-hand matrix */
	double *B = malloc(sizeof(double) * m * nrhs); /* Right-hand matrix */
	// double *X = malloc(sizeof(double) * n * nrhs);

	int i, base;
	    
	for (i = 0; i < num_pts; i++) {
	    base = 3 * i * n;

	    /* Row 1 */
	    A[base + 0] = (double) Vx(l_pts[i]);
	    A[base + 1] = (double) Vy(l_pts[i]);
	    A[base + 2] = (double) Vz(l_pts[i]);
	    A[base + 3] = 1.0;
	    A[base + 4] = A[base + 5] = A[base + 6] = A[base + 7] = 0.0;
	    A[base + 8] = A[base + 9] = A[base + 10] = A[base + 11] = 0.0;

	    base = (3 * i + 1) * n;

	    /* Row 2 */
	    A[base + 0] = A[base + 1] = A[base + 2] = A[base + 3] = 0.0;
	    A[base + 4] = (double) Vx(l_pts[i]);
	    A[base + 5] = (double) Vy(l_pts[i]);
	    A[base + 6] = (double) Vz(l_pts[i]);
	    A[base + 7] = 1.0;
	    A[base + 8] = A[base + 9] = A[base + 10] = A[base + 11] = 0.0;

	    base = (3 * i + 2) * n;

	    /* Row 2 */
	    A[base + 0] = A[base + 1] = A[base + 2] = A[base + 3] = 0.0;
	    A[base + 4] = A[base + 5] = A[base + 6] = A[base + 7] = 0.0;
	    A[base + 8] = (double) Vx(l_pts[i]);
	    A[base + 9] = (double) Vy(l_pts[i]);
	    A[base + 10] = (double) Vz(l_pts[i]);
	    A[base + 11] = 1.0;

	    B[3 * i + 0] = (double) Vx(r_pts[i]);
	    B[3 * i + 1] = (double) Vy(r_pts[i]);
	    B[3 * i + 2] = (double) Vz(r_pts[i]);
	}

	/* Run the driver to dgelsy */
	dgelsy_driver(A, B, Tout, m, n, nrhs);
	Tout[12] = 0.0; Tout[13] = 0.0; Tout[14] = 0.0; Tout[15] = 1.0;

	free(A);
	free(B);
	// free(X);
    }
}
