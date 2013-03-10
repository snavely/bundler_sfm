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

/* tps.c */
/* Routines for computing thin plate splines */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"
#include "tps.h"
#include "vector.h"

static double Ufn(double rsq) 
{
    if (rsq == 0.0)
	return 0.0;
    else
	return rsq * log(rsq);
}

/* Computes the TPS transformation that, when applied to the points
 * in l_pts, minimizes the least-squares error between the result and
 * the corresponding points in r_pts.
 * 
 * n -- number of points
 * r_pts -- matches
 * l_pts -- initial points 
 * Tout -- on return, contains the 3x3 transformation matrix */
void align_tps(int num_pts, v3_t *r_pts, v3_t *l_pts,double lambda, 
	       double *x_affine, double *x_weights,
	       double *y_affine, double *y_weights) 
{
    int num_eqns = num_pts + 3;
    int num_vars = num_pts + 3;
    
    double *A = (double *) malloc(sizeof(double) * num_eqns * num_vars);
    double *X = (double *) malloc(sizeof(double) * num_vars);
    double *B = (double *) malloc(sizeof(double) * num_eqns);

    int i, j;

    /* Fill in the top left part of A */
    for (i = 0; i < num_pts; i++) {
	for (j = 0; j < num_pts; j++) {
	    double dx = Vx(l_pts[i]) - Vx(l_pts[j]);
	    double dy = Vy(l_pts[i]) - Vy(l_pts[j]);

	    A[i * num_vars + j] = Ufn(dx * dx + dy * dy);

	    /* Add in the regularization */
	    if (i == j)
		A[i * num_vars + j] += lambda;
	}
    }

    /* Fill in the rest of A */
    for (i = 0; i < num_pts; i++) {
	A[i * num_vars + num_pts + 0] = 1.0;
	A[i * num_vars + num_pts + 1] = Vx(l_pts[i]);
	A[i * num_vars + num_pts + 2] = Vy(l_pts[i]);

	A[(num_pts + 0) * num_vars + i] = 1.0;
	A[(num_pts + 1) * num_vars + i] = Vx(l_pts[i]);
	A[(num_pts + 2) * num_vars + i] = Vy(l_pts[i]);
    }

    /* Zeroes elsewhere */
    for (i = num_pts; i < num_pts + 3; i++) {
	for (j = num_pts; j < num_pts + 3; j++) {
	    A[i * num_vars + j] = 0.0;
	}
    }

    /* Fill in B for the x parameters */
    for (i = 0; i < num_pts; i++) {
	B[i] = Vx(r_pts[i]);
    }

    for (i = num_pts; i < num_pts + 3; i++) {
	B[i]  = 0.0;
    }
    
    /* Solve */
    // dgelsy_driver(A, B, X, num_eqns, num_vars, 1);
    dgesv_driver(num_eqns, A, B, X);

    /* Copy out solution */
    memcpy(x_weights, X + 0, sizeof(double) * num_pts);
    x_affine[0] = X[num_pts + 1];
    x_affine[1] = X[num_pts + 2];
    x_affine[2] = X[num_pts + 0];



    /* Fill in B for the y parameters */
    for (i = 0; i < num_pts; i++) {
	B[i] = Vy(r_pts[i]);
    }

    for (i = num_pts; i < num_pts + 3; i++) {
	B[i]  = 0.0;
    }

    /* Solve */
    dgesv_driver(num_eqns, A, B, X);
    
    /* Copy out solution */
    memcpy(y_weights, X + 0, sizeof(double) * num_pts);
    y_affine[0] = X[num_pts + 1];
    y_affine[1] = X[num_pts + 2];
    y_affine[2] = X[num_pts + 0];

    free(A);
    free(B);
    free(X);
}

void align_tps_factored(int num_pts, v3_t *r_pts,
                        double *LU, int *ipiv, 
                        double *x_affine, double *x_weights,
                        double *y_affine, double *y_weights) 
{
    int num_eqns = num_pts + 3;
    int num_vars = num_pts + 3;

    double *X = (double *) malloc(sizeof(double) * num_vars);
    double *B = (double *) malloc(sizeof(double) * num_eqns);

    int i;

    /* Fill in B for the x parameters */
    for (i = 0; i < num_pts; i++) {
	B[i] = Vx(r_pts[i]);
    }

    for (i = num_pts; i < num_pts + 3; i++) {
	B[i]  = 0.0;
    }
    
    /* Solve */
    matrix_solve_lu(num_eqns, LU, ipiv, B, X);
    // dgesv_driver(num_eqns, A, B, X);

    /* Copy out solution */
    memcpy(x_weights, X + 0, sizeof(double) * num_pts);
    x_affine[0] = X[num_pts + 1];
    x_affine[1] = X[num_pts + 2];
    x_affine[2] = X[num_pts + 0];


    /* Fill in B for the y parameters */
    for (i = 0; i < num_pts; i++) {
	B[i] = Vy(r_pts[i]);
    }

    for (i = num_pts; i < num_pts + 3; i++) {
	B[i]  = 0.0;
    }

    /* Solve */
    matrix_solve_lu(num_eqns, LU, ipiv, B, X);
    
    /* Copy out solution */
    memcpy(y_weights, X + 0, sizeof(double) * num_pts);
    y_affine[0] = X[num_pts + 1];
    y_affine[1] = X[num_pts + 2];
    y_affine[2] = X[num_pts + 0];

    free(B);
    free(X);
}

void get_tps_basis(int num_pts, v3_t *l_pts, double lambda, 
                   double *LU, int *ipiv) 
{
    int num_eqns = num_pts + 3;
    int num_vars = num_pts + 3;
    
    double *A = (double *) malloc(sizeof(double) * num_eqns * num_vars);

    int i, j;

    /* Fill in the top left part of A */
    for (i = 0; i < num_pts; i++) {
	for (j = 0; j < num_pts; j++) {
	    double dx = Vx(l_pts[i]) - Vx(l_pts[j]);
	    double dy = Vy(l_pts[i]) - Vy(l_pts[j]);

	    A[i * num_vars + j] = Ufn(dx * dx + dy * dy);

	    /* Add in the regularization */
	    if (i == j)
		A[i * num_vars + j] += lambda;
	}
    }

    /* Fill in the rest of A */
    for (i = 0; i < num_pts; i++) {
	A[i * num_vars + num_pts + 0] = 1.0;
	A[i * num_vars + num_pts + 1] = Vx(l_pts[i]);
	A[i * num_vars + num_pts + 2] = Vy(l_pts[i]);

	A[(num_pts + 0) * num_vars + i] = 1.0;
	A[(num_pts + 1) * num_vars + i] = Vx(l_pts[i]);
	A[(num_pts + 2) * num_vars + i] = Vy(l_pts[i]);
    }

    /* Zeroes elsewhere */
    for (i = num_pts; i < num_pts + 3; i++) {
	for (j = num_pts; j < num_pts + 3; j++) {
	    A[i * num_vars + j] = 0.0;
	}
    }

    /* Factor the matrix */
    matrix_lu_no_transpose(num_eqns, A, LU, ipiv);

    free(A);
}
