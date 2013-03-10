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

/* horn.c */
/* Compute the closed-form least-squares solution to a rigid body alignment */

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "horn.h"
#include "matrix.h"
#include "qsort.h"
#include "svd.h"
#include "vector.h"

#ifdef WIN32
#define isnan _isnan
//#define isinf _isinf
#endif

/* Computes the closed-form least-squares solution to a rigid
 * body alignment.
 *
 * n: the number of points
 * right_pts: Target set of n points 
 * left_pts:  Source set of n points */
double align_horn(int n, v3_t *right_pts, v3_t *left_pts, 
		  double *R, double *T, 
		  double *Tout, double *scale, double *weight) {
    int i;
    v3_t right_centroid = v3_new(0.0, 0.0, 0.0);
    v3_t left_centroid = v3_new(0.0, 0.0, 0.0);
    double M[2][2] = { { 0.0, 0.0 }, 
                       { 0.0, 0.0 } };
    double MT[2][2];
    double MTM[2][2];
    double eval[2], sqrteval[2];
    double evec[2][2];
    double S[2][2], Sinv[2][2], U[2][2];
    double Tcenter[3][3] = { { 1.0, 0.0, 0.0 },
			     { 0.0, 1.0, 0.0 },
			     { 0.0, 0.0, 1.0 } };

    double Ttmp[3][3];

    double sum_num, sum_den, RMS_sum;

#if 1
    double weight_sum = 0.0;

    if (weight == NULL) {
        weight_sum = n;

        for (i = 0; i < n; i++) {
            right_centroid = 
                v3_add(right_centroid, right_pts[i]);
            left_centroid = 
                v3_add(left_centroid, left_pts[i]);
        }
        
        right_centroid = v3_scale(1.0 / weight_sum, right_centroid);
        left_centroid = v3_scale(1.0 / weight_sum, left_centroid);        
    } else {
        /* Compute the weighted centroid of both point sets */
        for (i = 0; i < n; i++) {
            right_centroid = 
                v3_add(right_centroid, v3_scale(weight[i], right_pts[i]));
            left_centroid = 
                v3_add(left_centroid, v3_scale(weight[i], left_pts[i]));
            weight_sum += weight[i];
        }

        right_centroid = v3_scale(1.0 / weight_sum, right_centroid);
        left_centroid = v3_scale(1.0 / weight_sum, left_centroid);
    }
#else
    /* Calculate the centroid of both sets of points */
    for (i = 0; i < n; i++) {
        right_centroid = v3_add(right_centroid, right_pts[i]);
        left_centroid = v3_add(left_centroid, left_pts[i]);
    }

    right_centroid = v3_scale(1.0 / n, right_centroid);
    left_centroid = v3_scale(1.0 / n, left_centroid);
#endif

    /* Compute the scale */
    sum_num = sum_den = 0.0;

    for (i = 0; i < n; i++) {
        v3_t r = v3_sub(right_centroid, right_pts[i]);
        v3_t l = v3_sub(left_centroid, left_pts[i]);
	
	sum_num = v3_magsq(r);
	sum_den = v3_magsq(l);
    }

    *scale = sqrt(sum_num / sum_den);

    /* Fill in the matrix M */
    for (i = 0; i < n; i++) {
        v3_t r = v3_sub(right_centroid, right_pts[i]);
        v3_t l = v3_sub(left_centroid, left_pts[i]);

	if (weight != NULL) {
	    M[0][0] += Vx(r) * Vx(l);
	    M[0][1] += Vx(r) * Vy(l);
	    M[1][0] += Vy(r) * Vx(l);
	    M[1][1] += Vy(r) * Vy(l);
	} else {
	    M[0][0] += Vx(r) * Vx(l);
	    M[0][1] += Vx(r) * Vy(l);
	    M[1][0] += Vy(r) * Vx(l);
	    M[1][1] += Vy(r) * Vy(l);
	}
    }

    /* Compute MTM */
    matrix_transpose(2, 2, (double *)M, (double *)MT);
    matrix_product(2, 2, 2, 2, (double *)MT, (double *)M, (double *)MTM);

    /* Calculate Sinv, the inverse of the square root of MTM */
    dgeev_driver(2, (double *)MTM, (double *)evec, eval);
    
    /* MTM = eval[0] * evec[0]T * evec[0] + eval[1] * evec[1]T * evec[1] */
    /* S = sqrt(eval[0]) * evec[0]T * evec[0] + sqrt(eval[1]) * evec[1]T * evec[1] */
    sqrteval[0] = sqrt(eval[0]);
    sqrteval[1] = sqrt(eval[1]);

    S[0][0] = 
        (sqrteval[0]) * evec[0][0] * evec[0][0] +
        (sqrteval[1]) * evec[1][0] * evec[1][0];
    S[0][1] = 
        (sqrteval[0]) * evec[0][0] * evec[0][1] +
        (sqrteval[1]) * evec[1][0] * evec[1][1];
    S[1][0] = 
        (sqrteval[0]) * evec[0][1] * evec[0][0] +
        (sqrteval[1]) * evec[1][1] * evec[1][0];
    S[1][1] = 
        (sqrteval[0]) * evec[0][1] * evec[0][1] +
        (sqrteval[1]) * evec[1][1] * evec[1][1];
    
    Sinv[0][0] = 
        (1.0 / sqrteval[0]) * evec[0][0] * evec[0][0] +
        (1.0 / sqrteval[1]) * evec[1][0] * evec[1][0];
    Sinv[0][1] = 
        (1.0 / sqrteval[0]) * evec[0][0] * evec[0][1] +
        (1.0 / sqrteval[1]) * evec[1][0] * evec[1][1];
    Sinv[1][0] = 
        (1.0 / sqrteval[0]) * evec[0][1] * evec[0][0] +
        (1.0 / sqrteval[1]) * evec[1][1] * evec[1][0];
    Sinv[1][1] = 
        (1.0 / sqrteval[0]) * evec[0][1] * evec[0][1] +
        (1.0 / sqrteval[1]) * evec[1][1] * evec[1][1];
    
    // matrix_product(2, 2, 2, 2, (double *)S, (double *)Sinv, (double *)U);

    /* U = M * Sinv */
    matrix_product(2, 2, 2, 2, (double *)M, (double *)Sinv, (double *)U);

    /* Fill in the rotation matrix */
    R[0] = U[0][0]; R[1] = U[0][1]; R[2] = 0.0;
    R[3] = U[1][0], R[4] = U[1][1]; R[5] = 0.0;
    R[6] = 0.0;     R[7] = 0.0;     R[8] = 1.0;

    // memcpy(R, U, sizeof(double) * 4);

    /* Fill in the translation matrix */
    T[0] = T[4] = T[8] = 1.0;
    T[1] = T[3] = T[6] = T[7] = 0.0;
    T[2] = Vx(right_centroid);
    T[5] = Vy(right_centroid);

    Tcenter[0][0] = *scale;
    Tcenter[1][1] = *scale;
    Tcenter[0][2] = -*scale * Vx(left_centroid);
    Tcenter[1][2] = -*scale * Vy(left_centroid);

    matrix_product(3, 3, 3, 3, T, R, (double *)Ttmp);

#if 0
#if 0
    /* Do the scaling */
    Ttmp[0][0] *= *scale;
    Ttmp[0][1] *= *scale;
    Ttmp[0][2] *= *scale;
    Ttmp[1][0] *= *scale;
    Ttmp[1][1] *= *scale;
    Ttmp[1][2] *= *scale;
#else
    Tcenter[0][0] *= *scale;
    Tcenter[0][2] *= *scale;
    Tcenter[1][1] *= *scale;
    Tcenter[1][2] *= *scale;
#endif
#endif

    matrix_product(3, 3, 3, 3, (double *)Ttmp, (double *)Tcenter, Tout);

    T[2] = Vx(v3_sub(right_centroid, left_centroid));
    T[5] = Vy(v3_sub(right_centroid, left_centroid));


    /* Now compute the RMS error between the points */
    RMS_sum = 0.0;

    for (i = 0; i < n; i++) {
        v3_t r = v3_sub(right_centroid, right_pts[i]);
        v3_t l = v3_sub(left_centroid, left_pts[i]);
	v3_t resid;

	/* Rotate, scale l */
	v3_t Rl, SRl;

	Vx(Rl) = R[0] * Vx(l) + R[1] * Vy(l) + R[2] * Vz(l);
	Vy(Rl) = R[3] * Vx(l) + R[4] * Vy(l) + R[5] * Vz(l);
	Vz(Rl) = R[6] * Vx(l) + R[7] * Vy(l) + R[8] * Vz(l);

	SRl = v3_scale(*scale, Rl);
	
	resid = v3_sub(r, SRl);
	RMS_sum += v3_magsq(resid);
    }
    
    return sqrt(RMS_sum / n);
}

/* Computes the closed-form least-squares solution to a rigid
 * body alignment.
 *
 * n: the number of points
 * right_pts: Target set of n points 
 * left_pts:  Source set of n points */
double align_horn_3D(int n, v3_t *right_pts, v3_t *left_pts, int scale_xform, 
		     double *Tout) {
    int i;
    v3_t right_centroid = v3_new(0.0, 0.0, 0.0);
    v3_t left_centroid = v3_new(0.0, 0.0, 0.0);
    double M[3][3] = { { 0.0, 0.0, 0.0, }, 
                       { 0.0, 0.0, 0.0, },
		       { 0.0, 0.0, 0.0, } };
    double MT[3][3];
    double MTM[3][3];
    double eval[3], sqrteval_inv[3];
    double evec[3][3], evec_tmp[3][3];
    double Sinv[3][3], U[3][3];
    double Tcenter[4][4] = { { 1.0, 0.0, 0.0, 0.0 },
			     { 0.0, 1.0, 0.0, 0.0 },
			     { 0.0, 0.0, 1.0, 0.0 },
			     { 0.0, 0.0, 0.0, 1.0 } };

    double Ttmp[4][4];
    double T[16], R[16];

    double sum_num, sum_den, scale, RMS_sum;

    int perm[3];

    /* Compute the centroid of both point sets */
    right_centroid = v3_mean(n, right_pts);
    left_centroid  = v3_mean(n, left_pts);

    /* Compute the scale */
    sum_num = sum_den = 0.0;

    for (i = 0; i < n; i++) {
        v3_t r = v3_sub(right_centroid, right_pts[i]);
        v3_t l = v3_sub(left_centroid, left_pts[i]);
	
	sum_num += v3_magsq(r);
	sum_den += v3_magsq(l);
    }

    scale = sqrt(sum_num / sum_den);

    /* Fill in the matrix M */
    for (i = 0; i < n; i++) {
        v3_t r = v3_sub(right_centroid, right_pts[i]);
        v3_t l = v3_sub(left_centroid, left_pts[i]);

	M[0][0] += Vx(r) * Vx(l);
	M[0][1] += Vx(r) * Vy(l);
	M[0][2] += Vx(r) * Vz(l);

	M[1][0] += Vy(r) * Vx(l);
	M[1][1] += Vy(r) * Vy(l);
	M[1][2] += Vy(r) * Vz(l);

	M[2][0] += Vz(r) * Vx(l);
	M[2][1] += Vz(r) * Vy(l);
	M[2][2] += Vz(r) * Vz(l);
    }

    /* Compute MTM */
    matrix_transpose(3, 3, (double *)M, (double *)MT);
    matrix_product(3, 3, 3, 3, (double *)MT, (double *)M, (double *)MTM);

    /* Calculate Sinv, the inverse of the square root of MTM */
    dgeev_driver(3, (double *)MTM, (double *)evec, eval);

    /* Sort the eigenvalues */
    qsort_descending();
    qsort_perm(3, eval, perm);
    
    memcpy(evec_tmp[0], evec[perm[0]], sizeof(double) * 3);
    memcpy(evec_tmp[1], evec[perm[1]], sizeof(double) * 3);
    memcpy(evec_tmp[2], evec[perm[2]], sizeof(double) * 3);
    memcpy(evec, evec_tmp, sizeof(double) * 9);

    sqrteval_inv[0] = 1.0 / sqrt(eval[0]);
    sqrteval_inv[1] = 1.0 / sqrt(eval[1]);

    if (eval[2] < 1.0e-8 * eval[0]) {
        sqrteval_inv[2] = 0.0;
    } else {
        sqrteval_inv[2] = 1.0 / sqrt(eval[2]);
    }

    Sinv[0][0] = 
        sqrteval_inv[0] * evec[0][0] * evec[0][0] +
        sqrteval_inv[1] * evec[1][0] * evec[1][0] + 
        sqrteval_inv[2] * evec[2][0] * evec[2][0];
    Sinv[0][1] = 
        sqrteval_inv[0] * evec[0][0] * evec[0][1] +
        sqrteval_inv[1] * evec[1][0] * evec[1][1] + 
        sqrteval_inv[2] * evec[2][0] * evec[2][1];
    Sinv[0][2] = 
        sqrteval_inv[0] * evec[0][0] * evec[0][2] +
        sqrteval_inv[1] * evec[1][0] * evec[1][2] + 
        sqrteval_inv[2] * evec[2][0] * evec[2][2];

    Sinv[1][0] = 
        sqrteval_inv[0] * evec[0][1] * evec[0][0] +
        sqrteval_inv[1] * evec[1][1] * evec[1][0] +
        sqrteval_inv[2] * evec[2][1] * evec[2][0];
    Sinv[1][1] = 
        sqrteval_inv[0] * evec[0][1] * evec[0][1] +
        sqrteval_inv[1] * evec[1][1] * evec[1][1] +
        sqrteval_inv[2] * evec[2][1] * evec[2][1];
    Sinv[1][2] = 
        sqrteval_inv[0] * evec[0][1] * evec[0][2] +
        sqrteval_inv[1] * evec[1][1] * evec[1][2] +
        sqrteval_inv[2] * evec[2][1] * evec[2][2];
    
    Sinv[2][0] = 
        sqrteval_inv[0] * evec[0][2] * evec[0][0] +
        sqrteval_inv[1] * evec[1][2] * evec[1][0] +
        sqrteval_inv[2] * evec[2][2] * evec[2][0];
    Sinv[2][1] = 
        sqrteval_inv[0] * evec[0][2] * evec[0][1] +
        sqrteval_inv[1] * evec[1][2] * evec[1][1] +
        sqrteval_inv[2] * evec[2][2] * evec[2][1];
    Sinv[2][2] = 
        sqrteval_inv[0] * evec[0][2] * evec[0][2] +
        sqrteval_inv[1] * evec[1][2] * evec[1][2] +
        sqrteval_inv[2] * evec[2][2] * evec[2][2];

    /* U = M * Sinv */
    matrix_product(3, 3, 3, 3, (double *)M, (double *)Sinv, (double *)U);

    if (eval[2] < 1.0e-8 * eval[0]) {    
        double u3u3[9], Utmp[9];
        matrix_transpose_product2(3, 1, 3, 1, evec[2], evec[2], u3u3);

        matrix_sum(3, 3, 3, 3, (double *) U, u3u3, Utmp);
        
        if (matrix_determinant3(Utmp) < 0.0) {
            printf("[align_horn_3D] Recomputing matrix...\n");
            matrix_diff(3, 3, 3, 3, (double *) U, u3u3, Utmp);
        }

        memcpy(U, Utmp, 9 * sizeof(double));
    }
    
    /* Fill in the rotation matrix */
    R[0]  = U[0][0]; R[1]  = U[0][1]; R[2]  = U[0][2]; R[3]  = 0.0;
    R[4]  = U[1][0]; R[5]  = U[1][1]; R[6]  = U[1][2]; R[7]  = 0.0;
    R[8]  = U[2][0]; R[9]  = U[2][1]; R[10] = U[2][2]; R[11] = 0.0;
    R[12] = 0.0;     R[13] = 0.0;     R[14] = 0.0;     R[15] = 1.0;
    
    /* Fill in the translation matrix */
    matrix_ident(4, T);
    T[3]  = Vx(right_centroid);
    T[7]  = Vy(right_centroid);
    T[11] = Vz(right_centroid);

    if (scale_xform == 0)
	scale = 1.0;

    Tcenter[0][0] = scale;
    Tcenter[1][1] = scale;
    Tcenter[2][2] = scale;
    
    Tcenter[0][3] = -scale * Vx(left_centroid);
    Tcenter[1][3] = -scale * Vy(left_centroid);
    Tcenter[2][3] = -scale * Vz(left_centroid);

    matrix_product(4, 4, 4, 4, T, R, (double *) Ttmp);
    matrix_product(4, 4, 4, 4, (double *)Ttmp, (double *)Tcenter, Tout);

#if 0
    T[2] = Vx(v3_sub(right_centroid, left_centroid));
    T[5] = Vy(v3_sub(right_centroid, left_centroid));
    T[8] = Vz(v3_sub(right_centroid, left_centroid));
#endif

    /* Now compute the RMS error between the points */
    RMS_sum = 0.0;

    for (i = 0; i < n; i++) {
	double left[4] = { Vx(left_pts[i]), 
			   Vy(left_pts[i]), 
			   Vz(left_pts[i]), 1.0 };
	double left_prime[3];
	double dx, dy, dz;

	matrix_product(4, 4, 4, 1, Tout, left, left_prime);

	dx = left_prime[0] - Vx(right_pts[i]);
	dy = left_prime[1] - Vy(right_pts[i]);
	dz = left_prime[2] - Vz(right_pts[i]);

	RMS_sum += dx * dx + dy * dy + dz * dz;
#if 0
        v3_t r = v3_sub(right_centroid, right_pts[i]);
        v3_t l = v3_sub(left_centroid, left_pts[i]);
	v3_t resid;

	/* Rotate, scale l */
	v3_t Rl, SRl;

	Vx(Rl) = R[0] * Vx(l) + R[1] * Vy(l) + R[2] * Vz(l);
	Vy(Rl) = R[3] * Vx(l) + R[4] * Vy(l) + R[5] * Vz(l);
	Vz(Rl) = R[6] * Vx(l) + R[7] * Vy(l) + R[8] * Vz(l);

	SRl = v3_scale(scale, Rl);
	
	resid = v3_sub(r, SRl);
	RMS_sum += v3_magsq(resid);
#endif
    }
    
    return sqrt(RMS_sum / n);
}


/* Computes the closed-form least-squares solution to a rigid
 * body alignment.
 *
 * n: the number of points
 * right_pts: Target set of n points 
 * left_pts:  Source set of n points */
double align_horn_3D_2(int n, v3_t *right_pts, v3_t *left_pts, int scale_xform,
                       double *Tout) 
{
    int i;
    v3_t right_centroid = v3_new(0.0, 0.0, 0.0);
    v3_t left_centroid = v3_new(0.0, 0.0, 0.0);
    double Tcenter[16] = { 1.0, 0.0, 0.0, 0.0,
                           0.0, 1.0, 0.0, 0.0,
                           0.0, 0.0, 1.0, 0.0,
                           0.0, 0.0, 0.0, 1.0 };

    double Ttmp[4][4];
    double T[16], R[16], R3x3[9];

    double sum_num, sum_den, scale, RMS_sum;

    v3_t *left_pts_zm  = malloc(sizeof(v3_t) * n);
    v3_t *right_pts_zm = malloc(sizeof(v3_t) * n);

    double error = 0.0;

    /* Compute the centroid of both point sets */
    right_centroid = v3_mean(n, right_pts);
    left_centroid  = v3_mean(n, left_pts);

    /* Compute the scale */
    sum_num = sum_den = 0.0;

    for (i = 0; i < n; i++) {
        v3_t r = v3_sub(right_centroid, right_pts[i]);
        v3_t l = v3_sub(left_centroid, left_pts[i]);

	sum_num += v3_magsq(r);
	sum_den += v3_magsq(l);
    }

    scale = sqrt(sum_num / sum_den);

    for (i = 0; i < n; i++) {
        v3_t r = v3_sub(right_centroid, right_pts[i]);
        v3_t l = v3_sub(left_centroid, left_pts[i]);

        right_pts_zm[i] = r;
        left_pts_zm[i]  = v3_scale(scale, l);
    }

    /* Compute the rotation */
    error = align_3D_rotation(n, right_pts_zm, left_pts_zm, R3x3);
    // printf("error[%d]: %0.3f\n", n, error);
    // matrix_print(3, 3, R3x3);

    /* Fill in the rotation matrix */
    R[0]  = R3x3[0]; R[1]  = R3x3[1]; R[2]  = R3x3[2]; R[3]  = 0.0;
    R[4]  = R3x3[3]; R[5]  = R3x3[4]; R[6]  = R3x3[5]; R[7]  = 0.0;
    R[8]  = R3x3[6]; R[9]  = R3x3[7]; R[10] = R3x3[8]; R[11] = 0.0;
    R[12] = 0.0;     R[13] = 0.0;     R[14] = 0.0;     R[15] = 1.0;
    
    /* Fill in the translation matrix */
    // matrix_ident(4, T);
    T[0] = 1.0;  T[1] = 0.0;  T[2] = 0.0;  T[3]  = Vx(right_centroid);
    T[4] = 0.0;  T[5] = 1.0;  T[6] = 0.0;  T[7]  = Vy(right_centroid);
    T[8] = 0.0;  T[9] = 0.0;  T[10] = 1.0; T[11] = Vz(right_centroid);
    T[12] = 0.0; T[13] = 0.0; T[14] = 0.0; T[15] = 1.0;

    if (scale_xform == 0)
	scale = 1.0;

    Tcenter[0] = scale;
    Tcenter[5] = scale;
    Tcenter[10] = scale;
    
    Tcenter[3] = -scale * Vx(left_centroid);
    Tcenter[7] = -scale * Vy(left_centroid);
    Tcenter[11] = -scale * Vz(left_centroid);

    matrix_product44(T, R, (double *) Ttmp);
    matrix_product44((double *)Ttmp, (double *)Tcenter, Tout);

    /* Now compute the RMS error between the points */
    RMS_sum = 0.0;

    for (i = 0; i < n; i++) {
	double left[4] = { Vx(left_pts[i]), 
			   Vy(left_pts[i]), 
			   Vz(left_pts[i]), 1.0 };
	double left_prime[3];
	double dx, dy, dz;

	matrix_product441(Tout, left, left_prime);

	dx = left_prime[0] - Vx(right_pts[i]);
	dy = left_prime[1] - Vy(right_pts[i]);
	dz = left_prime[2] - Vz(right_pts[i]);

	RMS_sum += dx * dx + dy * dy + dz * dz;
    }

    free(left_pts_zm);
    free(right_pts_zm);
    
    return sqrt(RMS_sum / n);
}

/* Align two sets of points with a 3D rotation */
double align_3D_rotation(int n, v3_t *r_pts, v3_t *l_pts, double *R)
{
    double A[9];
    double U[9], S[3], V[9], VT[9], RT[9];
    int i;
    double error;

#if 0
    if (n > 3) {
        printf("A:\n");
        for (i = 0; i < n; i++) {
            printf("%0.6f %0.6f %0.6f\n", 
                   Vx(r_pts[i]), Vy(r_pts[i]), Vz(r_pts[i]));
        }
        printf("B:\n");
        for (i = 0; i < n; i++) {
            printf("%0.6f %0.6f %0.6f\n", 
                   Vx(l_pts[i]), Vy(l_pts[i]), Vz(l_pts[i]));                
        }        
    }
#endif

    for (i = 0; i < 9; i++)
	A[i] = 0.0;

    for (i = 0; i < n; i++) {
        double *a = l_pts[i].p, *b = r_pts[i].p;
	// matrix_product(3, 1, 1, 3, l_pts[i].p, r_pts[i].p, tensor);
        A[0] += a[0] * b[0];
        A[1] += a[0] * b[1];
        A[2] += a[0] * b[2];

        A[3] += a[1] * b[0];
        A[4] += a[1] * b[1];
        A[5] += a[1] * b[2];

        A[6] += a[2] * b[0];
        A[7] += a[2] * b[1];
        A[8] += a[2] * b[2];
    }

    // dgesvd_driver(3, 3, A, U, S, VT);
    // printf("svd:\n");
    // matrix_print(3, 3, A);
    svd(3, 3, 1, 1, 1.0e-12, 1.0e-12, A, S, U, V, VT);
    
    // printf("U:\n");
    // matrix_print(3, 3, U);
    // printf("VT:\n");
    // matrix_print(3, 3, VT);
    // printf("S:\n");
    // matrix_print(3, 3, S);

    matrix_product33(U, VT, RT);
    matrix_transpose(3, 3, RT, R);

    // printf("R:\n");
    // matrix_print(3, 3, R);

    if (matrix_determinant3(R) < 0.0) {
        /* We're dealing with a reflection */
        double tmp[9];
        double reflectZ[9] = { 1.0, 0.0,  0.0,
                               0.0, 1.0,  0.0, 
                               0.0, 0.0, -1.0 };

        matrix_product33(U, reflectZ, tmp);
        matrix_product33(tmp, VT, RT);
        matrix_transpose(3, 3, RT, R);
    }

    /* Compute error */
    error = 0.0;
    for (i = 0; i < n; i++) {
	double rot[3];
	double diff[3];
	double dist;
	matrix_product331(R, l_pts[i].p, rot);
	matrix_diff(3, 1, 3, 1, rot, r_pts[i].p, diff);
	dist = matrix_norm(3, 1, diff);

        // printf("d[%d] = %0.6f\n", i, dist);
	error += dist;
    }

    return error / n;
}

double align_2D(int n, v3_t *right_pts, v3_t *left_pts, 
                double *R, double *T, 
                double *Tout, double *scale, double *weight) 
{
    int i;
    v3_t right_centroid = v3_mean(n, right_pts);
    v3_t left_centroid = v3_mean(n, left_pts);
    double Ttmp[3][3], Tcenter[3][3];

    /* Setup the matrix */
    int num_vars = 2;
    int num_eqns = 2 * n;
    double *A = (double *) malloc(sizeof(double) * num_vars * num_eqns);
    double *b = (double *) malloc(sizeof(double) * num_eqns);
    double x[2];

    double RMS_sum = 0.0;

    assert(n >= 2);
    for (i = 0; i < n; i++) {
        double *row1 = A + 2 * i * num_vars;
        double *row2 = A + (2 * i + 1) * num_vars;
        
        row1[0] = Vx(left_pts[i]) - Vx(left_centroid);
        row1[1] = Vy(left_pts[i]) - Vy(left_centroid);
        b[2 * i + 0] = Vx(right_pts[i]) - Vx(right_centroid);

        row2[0] = Vy(left_pts[i]) - Vy(left_centroid);
        row2[1] = -(Vx(left_pts[i]) - Vx(left_centroid));
        b[2 * i + 1] = Vy(right_pts[i]) - Vy(right_centroid);
    }

    /* Solve the system */
    dgelsy_driver(A, b, x, num_eqns, num_vars, 1);    

    /* Fill in the rotation matrix */
    R[0] = x[0];  R[1] = x[1]; R[2] = 0.0;
    R[3] = -x[1]; R[4] = x[0]; R[5] = 0.0;
    R[6] = 0.0;   R[7] = 0.0;  R[8] = 1.0;

    /* Fill in the translation matrix */
    matrix_ident(3, T);
    T[2] = Vx(right_centroid);
    T[5] = Vy(right_centroid);

    matrix_ident(3, (double *) Tcenter);
    Tcenter[0][2] = -Vx(left_centroid);
    Tcenter[1][2] = -Vy(left_centroid);

    matrix_product(3, 3, 3, 3, T, R, (double *)Ttmp);
    matrix_product(3, 3, 3, 3, (double *)Ttmp, (double *)Tcenter, Tout);

    T[2] = Vx(v3_sub(right_centroid, left_centroid));
    T[5] = Vy(v3_sub(right_centroid, left_centroid));

    /* Now compute the RMS error between the points */
    RMS_sum = 0.0;

    for (i = 0; i < n; i++) {
        v3_t r = v3_sub(right_pts[i], right_centroid);
        v3_t l = v3_sub(left_pts[i], left_centroid);
	v3_t resid;

	/* Rotate, scale l */
	v3_t SRl;

	Vx(SRl) = R[0] * Vx(l) + R[1] * Vy(l) + R[2] * Vz(l);
	Vy(SRl) = R[3] * Vx(l) + R[4] * Vy(l) + R[5] * Vz(l);
	Vz(SRl) = R[6] * Vx(l) + R[7] * Vy(l) + R[8] * Vz(l);

	// SRl = v3_scale(*scale, Rl);
	
	resid = v3_sub(r, SRl);
	RMS_sum += v3_magsq(resid);
    }

    free(A);
    free(b);
    
    return sqrt(RMS_sum / n);
}

/* Align two sets of points with a 2D similarity transform */
int align_2D_ransac(int n, v3_t *r_pts, v3_t *l_pts, 
                    int num_ransac_rounds, double ransac_thresh,
                    double *Tret)
{
    int round;    
#define MIN_SUPPORT 2
    v3_t *l_inliers, *r_inliers;
    int num_inliers, max_inliers = 0;
    double Tbest[9], R[9], T[9], scale;
    int i;

    if (n < 3) {
	printf("[align_2D_ransac] Error: need at least 3 points!\n");
	return 0;
    }

    l_inliers = (v3_t *) malloc(sizeof(v3_t) * n);
    r_inliers = (v3_t *) malloc(sizeof(v3_t) * n);

    for (round = 0; round < num_ransac_rounds; round++) {
	int support[MIN_SUPPORT];
	int i, j;
    v3_t r_mean, l_mean, r0, l0;
	v3_t r_pts_small[MIN_SUPPORT], l_pts_small[MIN_SUPPORT];
	double Rtmp[9], T1tmp[9], T2tmp[9], tmp[9], Tout[9];
    double a, b;

    for (i = 0; i < MIN_SUPPORT; i++) {
        /* Select an index from 0 to n-1 */
        int idx, reselect;

        do {
            reselect = 0;
            idx = rand() % n;
            for (j = 0; j < i; j++) {
                if (support[j] == idx) {
                    reselect = 1;
                    break;
                }
            }
        } while (reselect);

        support[i] = idx;
        r_pts_small[i] = r_pts[idx];
        l_pts_small[i] = l_pts[idx];
    }

    r_mean = v3_scale(0.5, v3_add(r_pts_small[0], r_pts_small[1]));
    l_mean = v3_scale(0.5, v3_add(l_pts_small[0], l_pts_small[1]));

    r0 = v3_sub(r_pts_small[0], r_mean);
    // v3_t r1 = v3_sub(r_pts_small[1], mean);
    l0 = v3_sub(l_pts_small[0], l_mean);
    // v3_t l1 = v3_sub(l_pts_small[1], mean);

    a = (Vy(r0) + Vx(r0) * Vx(l0) / Vy(l0)) / 
        (Vy(l0) + Vx(l0) * Vx(l0) / Vy(l0));
    b = (Vx(r0) - a * Vx(l0)) / Vy(l0);

    Rtmp[0] = a;  Rtmp[1] = b;  Rtmp[2] = 0.0;
    Rtmp[3] = -b; Rtmp[4] = a;  Rtmp[5] = 0.0;
    Rtmp[6] = 0;  Rtmp[7] = 0;  Rtmp[8] = 1.0;

    matrix_ident(3, T1tmp);
    T1tmp[2] = -Vx(l_mean);
    T1tmp[5] = -Vy(l_mean);

    matrix_ident(3, T2tmp);
    T2tmp[2] = Vx(r_mean);
    T2tmp[5] = Vy(r_mean);        

    matrix_product(3, 3, 3, 3, Rtmp, T1tmp, tmp);
    matrix_product(3, 3, 3, 3, T2tmp, tmp, Tout);

	/* Count inliers */
	num_inliers = 0;
	for (i = 0; i < n; i++) {
	    double Tp[3];
	    double diff[3];
	    double dist;
	    matrix_product(3, 3, 3, 1, Tout, l_pts[i].p, Tp);
	    matrix_diff(3, 1, 3, 1, Tp, r_pts[i].p, diff);
	    dist = matrix_norm(3, 1, diff);
	    
	    if (dist < ransac_thresh) {
		num_inliers++;
	    }

	    if (num_inliers > max_inliers) {
		max_inliers = num_inliers;
		memcpy(Tbest, Tout, sizeof(double) * 9);
                // printf(" inliers_new: %d\n", num_inliers);
	    }
	}
    }

#if 0
    /* Reestimate using all inliers */
    num_inliers = 0;
    for (i = 0; i < n; i++) {
	double Tp[3];
	double diff[3];
	double dist;
	matrix_product(3, 3, 3, 1, Tbest, l_pts[i].p, Tp);
	matrix_diff(3, 1, 3, 1, Tp, r_pts[i].p, diff);
	dist = matrix_norm(3, 1, diff);
	    
	if (dist < ransac_thresh) {
	    r_inliers[num_inliers] = r_pts[i];
	    l_inliers[num_inliers] = l_pts[i];
	    num_inliers++;
	}
    }

    // align_horn(num_inliers, r_inliers, l_inliers, R, T, Tret, &scale, NULL);
    align_2D(num_inliers, r_inliers, l_inliers, R, T, Tret, &scale, NULL);
#else
    memcpy(Tret, Tbest, 9 * sizeof(double));
#endif

    free(r_inliers);
    free(l_inliers);

    return num_inliers;
#undef MIN_SUPPORT
}

/* Align two sets of points with a 2D similarity transform */
int align_horn_ransac(int n, v3_t *r_pts, v3_t *l_pts, 
                      int num_ransac_rounds, double ransac_thresh,
                      double *Tret)
{
    int round;    
#define MIN_SUPPORT 3
    v3_t *l_inliers, *r_inliers;
    int num_inliers, max_inliers = 0;
    double Tbest[9], R[9], T[9], scale;
    int i;

    if (n < 3) {
	printf("[align_horn_ransac] Error: need at least 3 points!\n");
	return 0;
    }

    l_inliers = (v3_t *) malloc(sizeof(v3_t) * n);
    r_inliers = (v3_t *) malloc(sizeof(v3_t) * n);

    for (round = 0; round < num_ransac_rounds; round++) {
	int support[MIN_SUPPORT];
	int i, j;
	v3_t r_pts_small[MIN_SUPPORT], l_pts_small[MIN_SUPPORT];
	double Rtmp[9], Ttmp[9], Tout[9], scale_tmp;
	
	for (i = 0; i < MIN_SUPPORT; i++) {
	    /* Select an index from 0 to n-1 */
	    int idx, reselect;
	    do {
		reselect = 0;
		idx = rand() % n;
		for (j = 0; j < i; j++) {
		    if (support[j] == idx) {
			reselect = 1;
			break;
		    }
		}
	    } while (reselect);

	    support[i] = idx;
	    r_pts_small[i] = r_pts[idx];
	    l_pts_small[i] = l_pts[idx];
	}

        align_horn(MIN_SUPPORT, r_pts_small, l_pts_small, Rtmp, Ttmp, 
                   Tout, &scale_tmp, NULL);

	/* Count inliers */
	num_inliers = 0;
	for (i = 0; i < n; i++) {
	    double Tp[3];
	    double diff[3];
	    double dist;
	    matrix_product(3, 3, 3, 1, Tout, l_pts[i].p, Tp);
	    matrix_diff(3, 1, 3, 1, Tp, r_pts[i].p, diff);
	    dist = matrix_norm(3, 1, diff);
	    
	    if (dist < ransac_thresh) {
		num_inliers++;
	    }

	    if (num_inliers > max_inliers) {
		max_inliers = num_inliers;
		memcpy(Tbest, Tout, sizeof(double) * 9);
                // printf(" inliers_new: %d\n", num_inliers);
	    }
	}
    }

#if 0
    /* Reestimate using all inliers */
    num_inliers = 0;
    for (i = 0; i < n; i++) {
	double Tp[3];
	double diff[3];
	double dist;
	matrix_product(3, 3, 3, 1, Tbest, l_pts[i].p, Tp);
	matrix_diff(3, 1, 3, 1, Tp, r_pts[i].p, diff);
	dist = matrix_norm(3, 1, diff);
	    
	if (dist < ransac_thresh) {
	    r_inliers[num_inliers] = r_pts[i];
	    l_inliers[num_inliers] = l_pts[i];
	    num_inliers++;
	}
    }

    // printf(" inliers: %d\n", num_inliers);

    align_horn(num_inliers, r_inliers, l_inliers, R, T, Tret, &scale, NULL);
#else
    memcpy(Tret, Tbest, 9 * sizeof(double));
#endif

    free(r_inliers);
    free(l_inliers);

    return num_inliers;
#undef MIN_SUPPORT
}

/* Align two sets of points with a 3D similarity transform */
int align_horn_3D_ransac(int n, v3_t *r_pts, v3_t *l_pts, 
                         int num_ransac_rounds, double ransac_thresh,
                         double *Tret)
{
    int round;    
#define MIN_SUPPORT 3
    v3_t *l_inliers, *r_inliers;
    double *Vp, *TVp;
    int num_inliers, max_inliers = 0;
    double Tbest[16], TbestT[16];
    int i;
    double *ptr;

    double ransac_threshsq = ransac_thresh * ransac_thresh;

    if (n < MIN_SUPPORT) {
	printf("[align_horn_3D_ransac] Error: need at least %d points!\n",
               MIN_SUPPORT);
	return 0;
    }

    l_inliers = (v3_t *) malloc(sizeof(v3_t) * n);
    r_inliers = (v3_t *) malloc(sizeof(v3_t) * n);

    Vp = (double *) malloc(sizeof(double) * 4 * n);
    TVp = (double *) malloc(sizeof(double) * 4 * n);

    for (i = 0; i < n; i++) {
        memcpy(Vp + 4 * i, l_pts[i].p, 3 * sizeof(double));
        Vp[4 * i + 3] = 1.0;
    }

    for (round = 0; round < num_ransac_rounds; round++) {
	int support[MIN_SUPPORT];
	int i, j;
	v3_t r_pts_small[MIN_SUPPORT], l_pts_small[MIN_SUPPORT];
	double Tout[16], ToutT[16];
        int nan = 0;
	
	for (i = 0; i < MIN_SUPPORT; i++) {
	    /* Select an index from 0 to n-1 */
	    int idx, reselect;
	    do {
		reselect = 0;
		idx = rand() % n;
		for (j = 0; j < i; j++) {
		    if (support[j] == idx) {
			reselect = 1;
			break;
		    }
		}
	    } while (reselect);

	    support[i] = idx;
	    r_pts_small[i] = r_pts[idx];
	    l_pts_small[i] = l_pts[idx];
	}

        align_horn_3D_2(MIN_SUPPORT, r_pts_small, l_pts_small, 1, Tout);

#if 1
        for (i = 0; i < 16; i++) {
            if (isnan(Tout[i]) || Tout[i] != Tout[i]) {
                nan = 1;
                break;
            }
        }
        
        if (nan == 1)
            continue;
#endif

	/* Count inliers */
	num_inliers = 0;

#if 0
	for (i = 0; i < n; i++) {
	    double Tp[4];
	    double diff[3];
	    double dist;

            double p[4] = { l_pts[i].p[0], l_pts[i].p[1], l_pts[i].p[2], 1.0 };

	    matrix_product(4, 4, 4, 1, Tout, p, Tp);
	    matrix_diff(3, 1, 3, 1, Tp, r_pts[i].p, diff);
	    dist = matrix_norm(3, 1, diff);
	    
	    if (dist < ransac_thresh) {
		num_inliers++;
	    }
	}
#else
        matrix_transpose(4, 4, Tout, ToutT);
        matrix_product_old(n, 4, 4, 4, Vp, ToutT, TVp);

        ptr = TVp;

        for (i = 0; i < n; i++) {
            // double diff[3], dist;
            double dx, dy, dz, dist;
            // matrix_diff(3, 1, 3, 1, TVp + 4 * i, r_pts[i].p, diff);
            dx = ptr[0] - r_pts[i].p[0];
            dy = ptr[1] - r_pts[i].p[1];
            dz = ptr[2] - r_pts[i].p[2];

	    dist = dx * dx + dy * dy + dz * dz; // matrix_normsq(3, 1, diff);
	    
	    if (dist < ransac_threshsq)
		num_inliers++;

            ptr += 4;
        }
#endif

        if (num_inliers > max_inliers) {
            max_inliers = num_inliers;
            memcpy(Tbest, Tout, sizeof(double) * 16);
        }
    }

    /* Reestimate using all inliers */
#if 0
    matrix_transpose(4, 4, Tbest, TbestT);
    matrix_product(n, 4, 4, 4, Vp, TbestT, TVp);

    num_inliers = 0;
    for (i = 0; i < n; i++) {
	// double Tp[4];
	double diff[3];
	double dist;

        matrix_diff(3, 1, 3, 1, TVp + 4 * i, r_pts[i].p, diff);

        // double p[4] = { l_pts[i].p[0], l_pts[i].p[1], l_pts[i].p[2], 1.0 };

	// matrix_product(4, 4, 4, 1, Tbest, p, Tp);
	// matrix_diff(3, 1, 3, 1, Tp, r_pts[i].p, diff);
	dist = matrix_normsq(3, 1, diff);
	    
	if (dist < ransac_threshsq) {
	    r_inliers[num_inliers] = r_pts[i];
	    l_inliers[num_inliers] = l_pts[i];
	    num_inliers++;
	}
    }

    align_horn_3D_2(num_inliers, r_inliers, l_inliers, 1, Tret);

    // memcpy(Tret, Tbest, 16 * sizeof(double));

    if (isnan(Tret[0]) || Tret[0] != Tret[0]) {
        printf("[align_horn_3D_ransac] nan at end [num_inliers: %d], "
               "restoring old matrix\n", num_inliers);
        memcpy(Tret, Tbest, sizeof(double) * 16);
    }
#else
    memcpy(Tret, Tbest, sizeof(double) * 16);
#endif

    free(r_inliers);
    free(l_inliers);
    free(Vp);
    free(TVp);

    return max_inliers;
#undef MIN_SUPPORT
}

/* Align two sets of points with a 3D rotation */
int align_3D_rotation_ransac(int n, v3_t *r_pts, v3_t *l_pts, 
			     int num_ransac_rounds, double ransac_thresh,
			     double *R)
{
    int round;    
    double error = 0.0;
#define MIN_SUPPORT 3
    // const int min_support = 3;
    v3_t *l_inliers, *r_inliers;
    int num_inliers, max_inliers = 0;
    double Rbest[9];
    int i;

    if (n < 3) {
	printf("[align_3D_rotation_ransac] Error: need at least 3 points!\n");
	return 0;
    }

    l_inliers = (v3_t *) malloc(sizeof(v3_t) * n);
    r_inliers = (v3_t *) malloc(sizeof(v3_t) * n);

    for (round = 0; round < num_ransac_rounds; round++) {
	int support[MIN_SUPPORT];
	int i, j;
	v3_t r_pts_small[MIN_SUPPORT], l_pts_small[MIN_SUPPORT];
	double Rtmp[9];
	
	for (i = 0; i < MIN_SUPPORT; i++) {
	    /* Select an index from 0 to n-1 */
	    int idx, reselect;
	    do {
		reselect = 0;
		idx = rand() % n;
		for (j = 0; j < i; j++) {
		    if (support[j] == idx) {
			reselect = 1;
			break;
		    }
		}
	    } while (reselect);

	    support[i] = idx;
	    r_pts_small[i] = r_pts[idx];
	    l_pts_small[i] = l_pts[idx];
	}

	align_3D_rotation(MIN_SUPPORT, r_pts_small, l_pts_small, Rtmp);

	/* Count inliers */
	num_inliers = 0;

	for (i = 0; i < n; i++) {
	    double rot[3];
	    double diff[3];
	    double dist;
	    matrix_product(3, 3, 3, 1, Rtmp, l_pts[i].p, rot);
	    matrix_diff(3, 1, 3, 1, rot, r_pts[i].p, diff);
	    dist = matrix_norm(3, 1, diff);
	    
	    if (dist < ransac_thresh) {
		num_inliers++;
	    }

	    if (num_inliers > max_inliers) {
		max_inliers = num_inliers;
		memcpy(Rbest, Rtmp, sizeof(double) * 9);
	    }
	}
    }

#if 0
    /* Reestimate using all inliers */
    num_inliers = 0;
    for (i = 0; i < n; i++) {
	double rot[3];
	double diff[3];
	double dist;
	matrix_product(3, 3, 3, 1, Rbest, l_pts[i].p, rot);
	matrix_diff(3, 1, 3, 1, rot, r_pts[i].p, diff);
	dist = matrix_norm(3, 1, diff);
	    
	if (dist < ransac_thresh) {
	    r_inliers[num_inliers] = r_pts[i];
	    l_inliers[num_inliers] = l_pts[i];
	    num_inliers++;
	}
    }

    error = align_3D_rotation(num_inliers, r_inliers, l_inliers, R);    

    printf("[align_3D_rotation] Error: %0.3f\n", error);

    free(r_inliers);
    free(l_inliers);

    return num_inliers;
#else
    memcpy(R, Rbest, 9 * sizeof(double));

    free(r_inliers);
    free(l_inliers);

    return max_inliers;
#endif

#undef MIN_SUPPORT
}
