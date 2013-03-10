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

/* fit.c */
/* Routines for fitting various curves to data */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fit.h"
#include "matrix.h"
#include "vector.h"

/* Fit a line to a set of 2D points */
void fit_2D_line(int num_pts, v2_t *pts, double *params) 
{
    double *A, *b;
    int i;
    
    A = (double *) malloc(sizeof(double) * num_pts * 2);
    b = (double *) malloc(sizeof(double) * num_pts);
    
    for (i = 0; i < num_pts; i++) {
	A[2 * i + 0] = Vx(pts[i]);
	A[2 * i + 1] = Vy(pts[i]);
	
	b[i] = -1.0;
    }

    dgelsy_driver(A, b, params, num_pts, 2, 1);

    free(A);
    free(b);
}

/* Fit a line to a set of 2D points */
double fit_2D_line_orthogonal_regression(int num_pts, v2_t *pts, 
					 double *params) 
{
    double *A;
    int i;
    v2_t *pts_zero_mean;
    v2_t normal, vec;
    double error = 0.0;

    /* Compute the mean of the points */
    v2_t mean = v2_new(0.0, 0.0);
    
    for (i = 0; i < num_pts; i++)
	mean = v2_add(mean, pts[i]);

    mean = v2_scale(1.0 / num_pts, mean);

    /* Create a zero-mean set of points */
    pts_zero_mean = (v2_t *) malloc(sizeof(v2_t) * num_pts);
    
    for (i = 0; i < num_pts; i++)
	pts_zero_mean[i] = v2_sub(pts[i], mean);

    A = (double *) malloc(sizeof(double) * num_pts * 2);
    
    for (i = 0; i < num_pts; i++) {
	A[2 * i + 0] = Vx(pts_zero_mean[i]);
	A[2 * i + 1] = Vy(pts_zero_mean[i]);
    }

    matrix_minimum_unit_norm_solution(num_pts, 2, A, params);

    /* Set the third parameter as the dot product of the mean with
     * the normal */
    normal = v2_new(params[0], params[1]);
    params[2] = -v2_dotp(mean, normal);

    if (params[2] > 0.0) {
	params[0] = -params[0];
	params[1] = -params[1];
	params[2] = -params[2];
    }

    normal = v2_new(params[0], params[1]);
    vec = v2_scale(params[2], normal);

    for (i = 0; i < num_pts; i++) {
	v2_t diff = v2_sub(pts[i], vec);
	double dot = v2_dotp(diff, normal);
	
	error += dot * dot;
    }

    free(A);
    free(pts_zero_mean);

    return error;
}

/* Compute the distance between a line and a point */
double line_point_distance(double *line, v2_t pt) 
{
    return fabs(line[0] * Vx(pt) + line[1] * Vy(pt) + line[2]);
}

/* Fit a line using orthogonal regression and RANSAC */
double fit_2D_line_ortreg_ransac(int num_pts, v2_t *pts, 
                                 int num_ransac_rounds, 
                                 double ransac_threshold, 
                                 int *num_inliers_out, double *params)
{
#define MIN_SAMPLE_SIZE 2

    int *inliers = (int *) malloc(sizeof(int) * num_pts);
    int *best_inliers = (int *) malloc(sizeof(int) * num_pts);
    v2_t *best_pts;
    int i, j, k;

	int max_inliers = 0;
	double best_error = 0.0;
	double best_params[4];

	int num_inliers = 0;

    if (num_pts < MIN_SAMPLE_SIZE) {
	printf("[fit_2D_line_ortreg_ransac] "
	       "Error: num_pts < min_sample_size\n");

	free(inliers);
	free(best_inliers);
	
	return -1.0;
    }

    for (i = 0; i < num_ransac_rounds; i++) {
	double params_tmp[3];
	
	int samples[MIN_SAMPLE_SIZE];
	v2_t pts_tmp[MIN_SAMPLE_SIZE];

	int num_inliers = 0;
	double error = 0.0;
	
	for (j = 0; j < MIN_SAMPLE_SIZE; j++) {
	    int sample;
	    int retry;

	    do {
		retry = 0;
		sample = rand() % num_pts;

		for (k = 0; k < j; k++) {
		    if (sample == samples[k]) {
			retry = 1;
			break;
		    }
		}
	    } while (retry == 1);
	    
	    samples[j] = sample;
	    pts_tmp[j] = pts[sample];
	}

	fit_2D_line_orthogonal_regression(MIN_SAMPLE_SIZE, 
                                          pts_tmp, params_tmp);


	/* Count the inliers */
	for (j = 0; j < num_pts; j++) {
	    double dist = line_point_distance(params_tmp, pts[j]);
	    
	    if (dist < ransac_threshold) {
		error += dist;
		inliers[num_inliers] = j;
		num_inliers++;
	    }
	}

	if (num_inliers > max_inliers) {
	    max_inliers = num_inliers;
	    best_error = error;
	    memcpy(best_params, params, sizeof(double) * 3);
	    memcpy(best_inliers, inliers, sizeof(int) * num_inliers);
	}
    }

    if (max_inliers < 2) {
	printf("[fit_2D_line_ortreg_ransac] "
	       "Error: couldn't find enough inliers!\n");
    }

    /* Re-estimate the plane using all the inliers */
    best_pts = (v2_t *) malloc(sizeof(v2_t) * max_inliers);
    for (i = 0; i < max_inliers; i++) {	
	best_pts[i] = pts[best_inliers[i]];
    }

    fit_2D_line_orthogonal_regression(max_inliers, best_pts, params);

    num_inliers = 0;
    for (j = 0; j < num_pts; j++) {
        double dist = line_point_distance(params, pts[j]);
	    
        if (dist < ransac_threshold) {
            num_inliers++;
        }
    }

    printf(" %d / %d\n", max_inliers, num_inliers);

    free(inliers);
    free(best_inliers);
    free(best_pts);

    *num_inliers_out = max_inliers;

    return best_error / max_inliers;

#undef MIN_SAMPLE_SIZE
}


/* Fit a plane to a set of 3D points */
void fit_3D_plane(int num_pts, v3_t *pts, double *params) 
{
    double *A, *b;
    int i;
    
    A = (double *) malloc(sizeof(double) * num_pts * 3);
    b = (double *) malloc(sizeof(double) * num_pts);
    
    for (i = 0; i < num_pts; i++) {
	A[3 * i + 0] = Vx(pts[i]);
	A[3 * i + 1] = Vy(pts[i]);
	A[3 * i + 2] = Vz(pts[i]);
	
	b[i] = -1.0;
    }

    dgelsy_driver(A, b, params, num_pts, 3, 1);

    free(A);
    free(b);
}

/* Fit a plane (which passes through the origin) to a set of 3D points */
void fit_3D_plane_through_origin(int num_pts, v3_t *pts, double *params)
{
    double *A;
    int i;
    
    A = (double *) malloc(sizeof(double) * num_pts * 3);
    
    for (i = 0; i < num_pts; i++) {
	A[3 * i + 0] = Vx(pts[i]);
	A[3 * i + 1] = Vy(pts[i]);
	A[3 * i + 2] = Vz(pts[i]);
    }

    // dgelsy_driver(A, b, params, num_pts, 3, 1);
    matrix_minimum_unit_norm_solution(num_pts, 3, A, params);

    free(A);
}

/* Fit a plane (which passes through the origin) to a set of 3D points */
void fit_3D_plane_through_origin2(int num_pts, v3_t *pts, double *params)
{
    double *A;
    double A2[9];
    int i;
    
    A = (double *) malloc(sizeof(double) * num_pts * 3);
    
    for (i = 0; i < num_pts; i++) {
	A[3 * i + 0] = Vx(pts[i]);
	A[3 * i + 1] = Vy(pts[i]);
	A[3 * i + 2] = Vz(pts[i]);
    }

    matrix_transpose_product(num_pts, 3, num_pts, 3, A, A, (double *) A2);

    // dgelsy_driver(A, b, params, num_pts, 3, 1);
    matrix_minimum_unit_norm_solution(3, 3, A2, params);

    free(A);
}

/* Fit a plane using orthogonal regression */
double fit_3D_plane_orthogonal_regression(int num_pts, v3_t *pts, 
					  double *params)
{
    /* First subtract out the mean of the points */
    v3_t mean = v3_new(0.0, 0.0, 0.0);
    v3_t *pts_zero_mean;
    v3_t normal;
    int i;
    double error = 0.0;
    v3_t vec;

    for (i = 0; i < num_pts; i++) {
	Vx(mean) += Vx(pts[i]);
	Vy(mean) += Vy(pts[i]);
	Vz(mean) += Vz(pts[i]);
    }
    
    mean = v3_scale(1.0 / num_pts, mean);

    pts_zero_mean = (v3_t *) malloc(sizeof(v3_t) * num_pts);

    for (i = 0; i < num_pts; i++)
	pts_zero_mean[i] = v3_sub(pts[i], mean);

    // fit_3D_plane_through_origin(num_pts, pts_zero_mean, params);
    fit_3D_plane_through_origin2(num_pts, pts_zero_mean, params);
    
    /* Set the fourth parameter as the dot product of the mean with
       the normal */
    normal = v3_new(params[0], params[1], params[2]);
    params[3] = -v3_dotp(mean, normal);

    if (params[3] > 0.0) {
	params[0] = -params[0];
	params[1] = -params[1];
	params[2] = -params[2];
	params[3] = -params[3];
    }

    normal = v3_new(params[0], params[1], params[2]);
    vec = v3_scale(-params[3], normal);

    for (i = 0; i < num_pts; i++) {
	v3_t diff = v3_sub(pts[i], vec);
	double dot = v3_dotp(diff, normal);
	
	error += dot * dot;
    }

    free(pts_zero_mean);

    return error;
}

/* Project a point onto a plane */
v3_t project_point_onto_plane(double *plane, v3_t pt)
{
    v3_t normal = v3_new(plane[0], plane[1], plane[2]);
    v3_t vec = v3_scale(-plane[3], normal);
    v3_t diff = v3_sub(pt, vec);

    double dot = v3_dotp(diff, normal);
    v3_t par = v3_scale(dot, normal);
    v3_t perp = v3_sub(diff, par);
    
    return v3_add(perp, vec);
}

/* Compute the distance between a plane and a point */
double plane_point_distance(double *plane, v3_t pt) 
{
    double dist = 
        Vx(pt) * plane[0] + Vy(pt) * plane[1] + Vz(pt) * plane[2] + plane[3];

    return fabs(dist);
}

/* Fit a plane using orthogonal regression and RANSAC */
double fit_3D_plane_ortreg_ransac(int num_pts, v3_t *pts, 
				  int num_ransac_rounds, 
				  double ransac_threshold, 
				  int *num_inliers_out, double *params)
{
#define MIN_SAMPLE_SIZE 3
	// const int min_sample_size = 3; /* Three points determine a plane */
    int *inliers = (int *) malloc(sizeof(int) * num_pts);
    int *best_inliers = (int *) malloc(sizeof(int) * num_pts);
    v3_t *best_pts;
    int i, j, k;

	int max_inliers = 0;
	double best_error = 0.0;
	double best_params[4];

	int num_inliers = 0;

    if (num_pts < MIN_SAMPLE_SIZE) {
	printf("[fit_3D_plane_ortreg_ransac] "
	       "Error: num_pts < min_sample_size\n");

	free(inliers);
	free(best_inliers);
	
	return -1.0;
    }

    for (i = 0; i < num_ransac_rounds; i++) {
	double params_tmp[4];
	
	int samples[MIN_SAMPLE_SIZE];
	v3_t pts_tmp[MIN_SAMPLE_SIZE];

	int num_inliers = 0;
	double error = 0.0;
	
	for (j = 0; j < MIN_SAMPLE_SIZE; j++) {
	    int sample;
	    int retry;

	    do {
		retry = 0;
		sample = rand() % num_pts;

		for (k = 0; k < j; k++) {
		    if (sample == samples[k]) {
			retry = 1;
			break;
		    }
		}
	    } while (retry == 1);
	    
	    samples[j] = sample;
	    pts_tmp[j] = pts[sample];
	}

	fit_3D_plane_orthogonal_regression(MIN_SAMPLE_SIZE, 
					   pts_tmp, params_tmp);


	/* Count the inliers */
	for (j = 0; j < num_pts; j++) {
	    double dist = plane_point_distance(params_tmp, pts[j]);
	    
	    if (dist < ransac_threshold) {
		error += dist;
		inliers[num_inliers] = j;
		num_inliers++;
	    }
	}

	if (num_inliers > max_inliers) {
	    max_inliers = num_inliers;
	    best_error = error;
	    memcpy(best_params, params, sizeof(double) * 4);
	    memcpy(best_inliers, inliers, sizeof(int) * num_inliers);
	}
    }

    if (max_inliers < 3) {
	printf("[fit_3D_plane_ortreg_ransac] "
	       "Error: couldn't find enough inliers!\n");
    }

    /* Re-estimate the plane using all the inliers */
    best_pts = (v3_t *) malloc(sizeof(v3_t) * max_inliers);
    for (i = 0; i < max_inliers; i++) {	
	best_pts[i] = pts[best_inliers[i]];
    }

    fit_3D_plane_orthogonal_regression(max_inliers, best_pts, params);

    num_inliers = 0;
    for (j = 0; j < num_pts; j++) {
        double dist = plane_point_distance(params, pts[j]);
	    
        if (dist < ransac_threshold) {
            num_inliers++;
        }
    }

    // printf(" %d / %d\n", max_inliers, num_inliers);

    free(inliers);
    free(best_inliers);
    free(best_pts);

    *num_inliers_out = max_inliers;

    return best_error / max_inliers;

#undef MIN_SAMPLE_SIZE
}
