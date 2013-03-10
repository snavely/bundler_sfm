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

/* homography.h */
/* Computes a homography */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "homography.h"
#include "matrix.h"
#include "vector.h"

static v3_t *condition_points(int num_points, v3_t *pts, double *T) {
    v3_t *pts_new = (v3_t *) malloc(sizeof(v3_t) * num_points);

    v3_t mean = v3_mean(num_points, pts);
    double total_dist = 0.0;
    double avg_dist;
    double factor;
    int i;

    for (i = 0; i < num_points; i++) {
	double dx = Vx(pts[i]) - Vx(mean);
	double dy = Vy(pts[i]) - Vy(mean);
	total_dist += sqrt(dx * dx + dy * dy);
    }

    avg_dist = total_dist / num_points;
    factor = sqrt(2.0) / avg_dist;

    for (i = 0; i < num_points; i++) {
	double x = factor * (Vx(pts[i]) - Vx(mean));
	double y = factor * (Vy(pts[i]) - Vy(mean));
	pts_new[i] = v3_new(x, y, 1.0);
    }

    T[0] = factor;  T[1] = 0.0;     T[2] = -factor * Vx(mean);
    T[3] = 0.0;     T[4] = factor;  T[5] = -factor * Vy(mean);
    T[6] = 0.0;      T[7] = 0.0;    T[8] = 1.0;

    return pts_new;
}

/* Computes the homography that, when applied to the points in l_pts,
 * minimizes the least-squares error between the result and the
 * corresponding points in r_pts.
 * 
 * n -- number of points
 * r_pts -- matches
 * l_pts -- initial points 
 * Tout -- on return, contains the 3x3 transformation matrix */
void align_homography(int num_pts, v3_t *r_pts, v3_t *l_pts, 
		      double *Tout, int refine) 
{
    int m = num_pts * 2;
    int n = 8;
    int nrhs = 1;
    int i, base;

    double *A = malloc(sizeof(double) * m * n);    /* Left-hand matrix */
    double *B = malloc(sizeof(double) * m * nrhs); /* Right-hand matrix */

    double Ttmp[9];
    double T1[9], T2[9];

#define _CONDITION_
#ifdef _CONDITION_
    /* Normalize the points */
    v3_t *r_pts_norm = condition_points(num_pts, r_pts, T1);
    v3_t *l_pts_norm = condition_points(num_pts, l_pts, T2);
    double T1inv[9];
#else
    v3_t *r_pts_norm = r_pts;
    v3_t *l_pts_norm = l_pts;
#endif

    for (i = 0; i < num_pts; i++) {
	base = 2 * i * n;
	A[base + 0] = Vx(l_pts_norm[i]);
	A[base + 1] = Vy(l_pts_norm[i]);
	A[base + 2] = 1.0;
	A[base + 3] = A[base + 4] = A[base + 5] = 0.0;
	A[base + 6] = -Vx(l_pts_norm[i]) * Vx(r_pts_norm[i]);
	A[base + 7] = -Vy(l_pts_norm[i]) * Vx(r_pts_norm[i]);

	base = (2 * i + 1) * n;
	A[base + 0] = A[base + 1] = A[base + 2] = 0.0;
	A[base + 3] = Vx(l_pts_norm[i]);
	A[base + 4] = Vy(l_pts_norm[i]);
	A[base + 5] = 1.0;
	A[base + 6] = -Vx(l_pts_norm[i]) * Vy(r_pts_norm[i]);
	A[base + 7] = -Vy(l_pts_norm[i]) * Vy(r_pts_norm[i]);
    
	B[2 * i + 0] = Vx(r_pts_norm[i]);
	B[2 * i + 1] = Vy(r_pts_norm[i]);
    }

    /* Make the call to dgelsy */
    dgelsy_driver(A, B, Tout, m, n, nrhs);

    Tout[8] = 1.0;

#ifdef _CONDITION_
    /* Undo normalization */
    matrix_invert(3, T1, T1inv);

    matrix_product(3, 3, 3, 3, T1inv, Tout, Ttmp);
    matrix_product(3, 3, 3, 3, Ttmp, T2, Tout);

    matrix_scale(3, 3, Tout, 1.0 / Tout[8], Tout);
#endif

    if (refine) {
	memcpy(Ttmp, Tout, sizeof(double) * 9);
	align_homography_non_linear(num_pts, r_pts, l_pts, Ttmp, Tout);
    }
    
    free(A);
    free(B);

#ifdef _CONDITION_
    free(r_pts_norm);
    free(l_pts_norm);
#endif
}

static int global_num_pts;
static v3_t *global_r_pts;
static v3_t *global_l_pts;
static int global_round;

static void homography_resids(int *m, int *n, double *x, double *fvec, int *iflag)
{
    int i;

    double resids = 0.0;

    double H[9];
    memcpy(H, x, 8 * sizeof(double));
    H[8] = 1.0;

    if (*iflag == 0 && global_num_pts > 4) {
	printf("[Round %d]\n", global_round);
	printf("  H=(%0.5f, %0.5f, %0.5f, %0.5f, %0.5f, %0.5f, %0.5f, %0.5f, %0.1f)\n", 
	       H[0], H[1], H[2], H[3], H[4], H[5], H[6], H[7], H[8]);
	global_round++;
    }
    
    for (i = 0; i < global_num_pts; i++) {
	double p[3], q[3];

	p[0] = Vx(global_l_pts[i]);
	p[1] = Vy(global_l_pts[i]);
	p[2] = Vz(global_l_pts[i]);

	if (*iflag == 0 && global_num_pts > 4)
	    printf("    p=(%0.3f, %0.3f, %0.3f)\n", p[0], p[1], p[2]);

	matrix_product(3, 3, 3, 1, H, p, q);

	if (*iflag == 0 && global_num_pts > 4)
	    printf("    q=(%0.3f, %0.3f, %0.3f)\n", q[0], q[1], q[2]);
	
	q[0] /= q[2];
	q[1] /= q[2];
	
	fvec[2 * i + 0] = q[0] - Vx(global_r_pts[i]);
	fvec[2 * i + 1] = q[1] - Vy(global_r_pts[i]);

	if (*iflag == 0 && global_num_pts > 4)
	    printf("    (%0.3f, %0.3f) ==> (%0.3f, %0.3f)\n", q[0], q[1], Vx(global_r_pts[i]), Vy(global_r_pts[i]));
    }

    for (i = 0; i < 2 * global_num_pts; i++) {
	resids += fvec[i] * fvec[i];
    }
    
    if (*iflag == 0 && global_num_pts > 4) 
	printf("resids = %0.3f\n", resids);
}

/* Use non-linear least squares to refine a homography */
void align_homography_non_linear(int num_pts, v3_t *r_pts, v3_t *l_pts, 
				 double *Tin, double *Tout) 
{
    double x[8];

    if (num_pts > 4) {
	printf("pre: ");
	matrix_print(3, 3, Tin);
    }
    
    memcpy(x, Tin, 8 * sizeof(double));

    global_num_pts = num_pts;
    global_r_pts = r_pts;
    global_l_pts = l_pts;
    global_round = 0;

    lmdif_driver(homography_resids, 2 * num_pts, 8, x, 1.0e-4);
    
    memcpy(Tout, x, 8 * sizeof(double));
    Tout[8] = 1.0;

    if (num_pts > 4) {
	printf("post: ");
	matrix_print(3, 3, Tout);
    }
}

#if 0
/* Use RANSAC to estimate a homography */
void align_homography_ransac_matches(int num_pts, v3_t *a_pts, v3_t *b_pts, 
				     int num_trials, double threshold, 
                                     double success_ratio,
                                     int essential, double *H) 
{
#define MIN_SAMPLES 4
    int i, j, k, idx;

    v3_t l_pts_best[MIN_SAMPLES], r_pts_best[MIN_SAMPLES];

    double Hbest[MIN_SAMPLES];
    double *resid;
    double error_min;
    int inliers_max;

    double *a_matrix, *b_matrix;

    // double threshold = 1.0e-10;

    // srand(time(0));

    /* Make an array of all good correspondences */
    if (num_pts < MIN_SAMPLES) {
	printf("[align_homography_ransac] "
               "Could not find 8 good correspondences, "
	       "homography estimation failed\n");
	return;
    }

    a_matrix = malloc(sizeof(double) * 3 * num_pts);
    b_matrix = malloc(sizeof(double) * 3 * num_pts);

    for (i = 0; i < num_pts; i++) {
        a_matrix[i] = Vx(a_pts[i]);
        a_matrix[i+num_pts] = Vy(a_pts[i]);
        a_matrix[i+2*num_pts] = Vz(a_pts[i]);

        b_matrix[i] = Vx(b_pts[i]);
        b_matrix[i+num_pts] = Vy(b_pts[i]);
        b_matrix[i+2*num_pts] = Vz(b_pts[i]);        
    }

    error_min = DBL_MAX;
    inliers_max = 0;
    resid = (double *) malloc(sizeof(double) * num_pts);

    /* Estimate the homography using RANSAC */
    for (i = 0; i < num_trials; i++) {
	int idxs[MIN_SAMPLES];
	v3_t l_pts[MIN_SAMPLES], r_pts[MIN_SAMPLES];
	double Htmp[9];
	// double error;
	int num_inliers = 0;
	int success, nan = 0;
        int round = 0;

	/* Sample 4 random correspondences */
	for (j = 0; j < MIN_SAMPLES; j++) {
	    int reselect = 0;

            if (round == 1000)
                return;

	    idx = rand() % num_pts;
	    
	    /* Make sure we didn't sample this index yet */
	    for (k = 0; k < j; k++) {
		if (idx == idxs[k] ||
                    (Vx(a_pts[idx]) == Vx(a_pts[idxs[k]]) &&
                     Vy(a_pts[idx]) == Vy(a_pts[idxs[k]]) &&
                     Vz(a_pts[idx]) == Vz(a_pts[idxs[k]])) ||
                    (Vx(b_pts[idx]) == Vx(b_pts[idxs[k]]) &&
                     Vy(b_pts[idx]) == Vy(b_pts[idxs[k]]) &&
                     Vz(b_pts[idx]) == Vz(b_pts[idxs[k]]))) {
		    reselect = 1;
		    break;
		}
	    }

	    if (reselect) {
                round++;
		j--;
		continue;
	    }

	    idxs[j] = idx;
	}

	/* Fill in the left and right points */
	for (j = 0; j < 8; j++) {
	    l_pts[j] = b_pts[idxs[j]];
	    r_pts[j] = a_pts[idxs[j]];
	}

	/* Estimate the F-matrix */
	success = estimate_fmatrix_linear(8, r_pts, l_pts, essential, 
                                          Ftmp, e1_tmp, e2_tmp);

        if (success == 0)
            nan = 1;

	for (j = 0; j < 9; j++) {
	    if (Ftmp[j] != Ftmp[j] /* isnan(Ftmp[j]) */) {
		printf("[estimate_fmatrix_ransac_matches] nan encountered\n");
                nan = 1;
		break;
	    }
	}

        /* Check for nan entries */
        if (isnan(Ftmp[0]) || isnan(Ftmp[1]) || isnan(Ftmp[2]) ||
            isnan(Ftmp[3]) || isnan(Ftmp[4]) || isnan(Ftmp[5]) ||
            isnan(Ftmp[6]) || isnan(Ftmp[7]) || isnan(Ftmp[8])) {
            printf("[estimate_fmatrix_ransac_matches] "
                   "nan matrix encountered\n");
            nan = 1;
        }

	if (nan) {
	    // error = DBL_MAX;
	    num_inliers = 0;
	} else {
            // printf("%0.3f\n", Ftmp[0]);

	    /* Compute residuals */
#if 1
	    for (j = 0; j < num_pts; j++) {
		resid[j] = fmatrix_compute_residual(Ftmp, a_pts[j], b_pts[j]);
		if (resid[j] < threshold)
		    num_inliers++;
	    }
#else
            fmatrix_compute_residuals(num_pts, Ftmp, a_matrix, b_matrix,
                                      resid);

            for (j = 0; j < num_pts; j++) {
		if (resid[j] < threshold)
		    num_inliers++;                
            }
#endif

#if 0
	    /* Find the median */
	    error = median(num_pts, resid);

	    if (error < error_min) {
		error_min = error;
		memcpy(Fbest, Ftmp, sizeof(double) * 9);
		memcpy(l_pts_best, l_pts, sizeof(v3_t) * 8);
		memcpy(r_pts_best, r_pts, sizeof(v3_t) * 8);
	    }
#else
	    if (num_inliers > inliers_max) {
		inliers_max = num_inliers;
		memcpy(Fbest, Ftmp, sizeof(double) * 9);
		memcpy(l_pts_best, l_pts, sizeof(v3_t) * 8);
		memcpy(r_pts_best, r_pts, sizeof(v3_t) * 8);
	    }
#endif
	}
	
#if 0
	if (error < threshold)
	    break;
#endif

        if ((double) num_inliers / num_pts > success_ratio)
            break;
    }

    // printf("Minimum error: %0.5e\n", error_min);
    // printf("Maximum inliers: %d\n", inliers_max);

    // matrix_print(3, 3, Fbest);

    free(resid);

    /* Copy out the F-matrix */
    memcpy(F, Fbest, sizeof(double) * 9);

    free(a_matrix);
    free(b_matrix);
}
#endif
