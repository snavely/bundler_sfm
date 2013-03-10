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

/* vector.c */
/* Routines for dealing with vectors */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "defines.h"
#include "matrix.h"
#include "util.h"
#include "vector.h"

/* ******************* 2D vectors ******************* */

iv2_t iv2_new(int16_t x, int16_t y) {
    iv2_t v = { { x, y } };
    return v;
}

iv2_t iv2_add(iv2_t u, iv2_t v) {
    return iv2_new(Vx(u) + Vx(v), Vy(u) + Vy(v));
}

v2_t v2_new(double x, double y) {
    v2_t v = { { x, y } };
    return v;
}

v2_t v2_add(v2_t u, v2_t v) {
    return v2_new(Vx(u) + Vx(v),
		  Vy(u) + Vy(v));
}

iv2_t iv2_sub(iv2_t u, iv2_t v) {
    return iv2_new(Vx(u) - Vx(v),
		   Vy(u) - Vy(v));
}

v2_t v2_sub(v2_t u, v2_t v) {
    return v2_new(Vx(u) - Vx(v),
		  Vy(u) - Vy(v));
}

v2_t v2_scale(double s, v2_t v) {
    return v2_new(s * Vx(v), s * Vy(v));
}

double v2_norm(v2_t v) {
    return Vx(v) * Vx(v) + Vy(v) * Vy(v);
}

/* Return a unit vector in the same direction as v, v != 0 */
v2_t v2_unit(const v2_t v) {
    double mag = v2_norm(v);

    if (mag == 0) 
        return v;
    
    return v2_scale(1.0 / mag, v);
}

/* Dot product */
double v2_dotp(v2_t u, v2_t v) {
    return Vx(u) * Vx(v) + Vy(u) * Vy(v);
}

double v2_angle(v2_t u, v2_t v) {
    double angle = atan2(Vy(v),Vx(v)) - atan2(Vy(u),Vx(u));

    /* Get in range [-pi,pi] */
    if (angle > M_PI)
        angle -= 2.0 * M_PI;
    if (angle < -M_PI)
        angle += 2.0 * M_PI;

    return angle;
}

/* Compute the pair-wise minimum / maximum of two vectors */
v2_t v2_minimum(v2_t u, v2_t v) {
    return v2_new(MIN(Vx(u), Vx(v)),
		  MIN(Vy(u), Vy(v)));
}

v2_t v2_maximum(v2_t u, v2_t v) {   
    return v2_new(MAX(Vx(u), Vx(v)),
		  MAX(Vy(u), Vy(v)));
}

/* Compute the mean of a set of vectors */
v2_t v2_mean(int n, v2_t *v) {
    int i;
    v2_t mean = v2_new(0.0, 0.0);
    
    for (i = 0; i < n; i++) {
	mean = v2_add(mean, v[i]);
    }

    return v2_scale(1.0 / n, mean);
}

/* Compute the covariance of a set of vectors */
void v2_covariance(int n, v2_t *v, v2_t mean, double *cov) {
    int i;
    
    cov[0] = cov[1] = cov[2] = cov[3] = 0.0;
    for (i = 0; i < n; i++) {
        v2_t vzm = v2_sub(v[i], mean);

        double xy = Vx(vzm) * Vy(vzm);
        cov[0] += Vx(vzm) * Vx(vzm);
        cov[1] += xy;
        cov[2] += xy;
        cov[3] += Vy(vzm) * Vy(vzm);
    }

    cov[0] /= n;
    cov[1] /= n;
    cov[2] /= n;
    cov[3] /= n;
}

/* Compute the centroid of an array of 2D vectors */
v2_t v2_compute_centroid(v2_t *pts, int num_pts) {
    int i;
    v2_t centroid = v2_new(0.0, 0.0);

    for (i = 0; i < num_pts; i++) {
	Vx(centroid) += Vx(pts[i]);
	Vy(centroid) += Vy(pts[i]);
    }

    return v2_scale(1.0 / ((double) num_pts), centroid);
}

iv2_t iv2_compute_centroid(iv2_t *pts, int num_pts) {
    int i;
    v2_t centroid = v2_new(0.0, 0.0);
    iv2_t i_centroid;

    for (i = 0; i < num_pts; i++) {
	Vx(centroid) += (double) Vx(pts[i]);
	Vy(centroid) += (double) Vy(pts[i]);
    }

#define ROUND(x) (((x) < 0.0) ? (int) ((x) - 0.5) : (int) ((x) + 0.5))
    centroid = v2_scale(1.0 / ((double) num_pts), centroid);
    i_centroid = iv2_new((int16_t) ROUND(Vx(centroid)), 
			 (int16_t) ROUND(Vy(centroid)));
    
    return i_centroid;
}

/* ******************* 3D vectors ******************* */

iv3_t iv3_new(int16_t x, int16_t y, int16_t z) {
    iv3_t v = { { x, y, z } };
    return v;
}

iv3_t iv3_add(iv3_t u, iv3_t v) {
    return iv3_new(Vx(u) + Vx(v), Vy(u) + Vy(v), Vz(u) + Vz(v));
}

v3_t v3_new(double x, double y, double z) {
    v3_t v = { { x, y, z } };
    return v;
}

v3_t v3_add(const v3_t u, const v3_t v) {
    return v3_new(Vx(u) + Vx(v), Vy(u) + Vy(v), Vz(u) + Vz(v));
}

v3_t v3_sub(const v3_t u, const v3_t v) {
    return v3_new(Vx(u) - Vx(v), Vy(u) - Vy(v), Vz(u) - Vz(v));
}

v3_t v3_scale(double c, const v3_t v) {
    return v3_new(c * Vx(v), c * Vy(v), c * Vz(v));
}

double v3_magsq(const v3_t v) {
    return Vx(v) * Vx(v) + Vy(v) * Vy(v) + Vz(v) * Vz(v);
}

double v3_mag(const v3_t v) {
    return sqrt(v3_magsq(v));
}

/* Compute coordinate-wise min, max of two vectors */
v3_t v3_min(const v3_t u, const v3_t v) {
    return v3_new(MIN(Vx(u), Vx(v)),
		  MIN(Vy(u), Vy(v)),
		  MIN(Vz(u), Vz(v)));
}

v3_t v3_max(const v3_t u, const v3_t v) {
    return v3_new(MAX(Vx(u), Vx(v)),
		  MAX(Vy(u), Vy(v)),
		  MAX(Vz(u), Vz(v)));
}

v3_t v3_unit(const v3_t v) {
    double mag = v3_mag(v);

    if (mag == 0) 
        return v;
    
    return v3_scale(1.0 / mag, v);
}

/* Assumes unit_normal is normalized */
v3_t v3_project(const v3_t v, const v3_t unit_normal) {
    double dot = v3_dotp(v, unit_normal);
    v3_t par = v3_scale(dot, unit_normal);
    v3_t perp = v3_sub(v, par);
    return perp;
}

/* Scale the vector so that the 3rd coordinate is 1 */
v3_t v3_homogenize(const v3_t v) {
    if (Vz(v) == 0.0) {
	return v3_new(DBL_MAX, DBL_MAX, 1.0);
    } else {
	return v3_new(Vx(v) / Vz(v), Vy(v) / Vz(v), 1.0);
    }
}

double v3_dotp(const v3_t u, const v3_t v) {
    return Vx(u) * Vx(v) + Vy(u) * Vy(v) + Vz(u) * Vz(v);
}

v3_t v3_cross(const v3_t u, const v3_t v) {
    return v3_new(Vy(u) * Vz(v) - Vz(u) * Vy(v),
                  Vz(u) * Vx(v) - Vx(u) * Vz(v),
                  Vx(u) * Vy(v) - Vy(u) * Vx(v));
}

/* Compute the mean of a set of vectors */
v3_t v3_mean(int n, const v3_t *v) {
    int i;
    v3_t mean = v3_new(0.0, 0.0, 0.0);
    
    for (i = 0; i < n; i++) {
	mean = v3_add(mean, v[i]);
    }

    return v3_scale(1.0 / n, mean);
}

/* Compute the "median" of a set of vectors */
v3_t v3_median(int n, const v3_t *v) {
    int i;
    v3_t median = v3_new(0.0, 0.0, 0.0);
    double min_dist = DBL_MAX;

    for (i = 0; i < n; i++) {
        double dist = 0.0;
        int j;

        for (j = 0; j < n; j++) {
            dist += v3_mag(v3_sub(v[j], v[i]));
        }
        
        if (dist < min_dist) {
            min_dist = dist;
            median = v[i];
        }
    }

    return median;
}

double v3_variance_zm(int n, const v3_t *v)
{
    int i;
    double sum = 0.0;
    for (i = 0; i < n; i++) {
        sum += v3_magsq(v[i]);
    }

    return sum / n;
}

/* Compute the covariance of a set of (zero-mean) vectors */
void v3_covariance_zm(int n, const v3_t *v, double *cov) {
    int i;
    
    for (i = 0; i < 9; i++)
        cov[i] = 0.0;

    for (i = 0; i < n; i++) {
        double tmp[9];
        matrix_product(3, 1, 1, 3, v[i].p, v[i].p, tmp);
        matrix_sum(3, 3, 3, 3, cov, tmp, cov);
    }

    matrix_scale(3, 3, cov, 1.0 / n, cov);
}


/* Compute the mean of a set of vectors */
void v3_svd(int n, const v3_t *v, double *U, double *S, double *VT) {
    double A[9] = { 0.0, 0.0, 0.0,
		    0.0, 0.0, 0.0,
		    0.0, 0.0, 0.0 };
    int i;

    for (i = 0; i < n; i++) {
	double tensor[9];
	matrix_product(3, 1, 1, 3, v[i].p, v[i].p, tensor);
	matrix_sum(3, 3, 3, 3, A, tensor, A);
    }
    
    dgesvd_driver(3, 3, A, U, S, VT);
}

/* Find the vector in v that is furthest from u */
v3_t v3_extremum(int n, const v3_t *v, const v3_t u) {
    int idx = v3_extremum_idx(n, v, u);
    return v[idx];
}

v3_t v3_extremum2(int n, const v3_t *a, const v3_t u, v3_t v) {
    int i;
    int max_idx = -1;
    double max_dist = 0.0;

    for (i = 0; i < n; i++) {
	v3_t diff1 = v3_sub(a[i], u);
	v3_t diff2 = v3_sub(a[i], v);
	double distsq1 = v3_magsq(diff1);
	double distsq2 = v3_magsq(diff2);
	if (MIN(distsq1, distsq2) > max_dist) {
	    max_idx = i;
	    max_dist = MIN(distsq1, distsq2);
	}
    }

    if (max_idx == -1) {
	printf("[v3_extremum2] Couldn't find extremum!\n");
	return v3_new(0.0, 0.0, 0.0);
    }

    return a[max_idx];
}

int v3_extremum_idx(int n, const v3_t *v, const v3_t u) {
    int i;
    int max_idx = -1;
    double max_dist = 0.0;

    for (i = 0; i < n; i++) {
	v3_t diff = v3_sub(v[i], u);
	double distsq = v3_magsq(diff);
	if (distsq > max_dist) {
	    max_idx = i;
	    max_dist = distsq;
	}
    }

    if (max_idx == -1) {
	printf("[v3_extremum] Couldn't find extremum!\n");
	return 0;
    }

    return max_idx;
}

void v3_print(const v3_t v) {
    printf("(%0.3f, %0.3f, %0.3f)\n", Vx(v), Vy(v), Vz(v));
}


/* ******************* ND vectors ******************* */

vec_t vec_new(int d) {
    vec_t v;

    v.d = d;
    v.p = malloc(sizeof(double) * d);

    return v;
}

vec_t vec_new_set(int d, double val) {
    vec_t v = vec_new(d);
    int i;

    for (i = 0; i < d; i++)
        Vn(v, i) = val;
    
    return v;
}

vec_t vec_add(vec_t u, vec_t v) {
    vec_t result = { 0, NULL };
    int i;

    if (u.d != v.d) {
        printf("Error: added vectors must have same dimension\n");
        return result;
    }

    result = vec_new(u.d);

    for (i = 0; i < u.d; i++)
        result.p[i] = u.p[i] + v.p[i];

    return result;
}

vec_t vec_sub(vec_t u, vec_t v) {
    vec_t result = { 0, NULL };
    int i;

    if (u.d != v.d) {
        printf("Error: subtracted vectors must have same dimension\n");
        return result;
    }
    
    result = vec_new(u.d);

    for (i = 0; i < u.d; i++)
        result.p[i] = u.p[i] - v.p[i];

    return result;
}

/* Scale the given vector with the given scalar (changing the
 * given vector) */
void vec_scale_inplace(double c, vec_t v) {
    int i;
    
    for (i = 0; i < v.d; i++)
	Vn(v, i) *= c;
}

/* Compute the norm (length squared) of the given vector */
double vec_norm(vec_t v) {
    int i;
    double n = 0.0;

    for (i = 0; i < v.d; i++) 
	n += v.p[i] * v.p[i];
    
    return n;
}

/* Copy one vector into another */
void vec_copy(vec_t dest, vec_t src) {
    int i;
    
    if (dest.d != src.d) {
	printf("[vec_copy] Mismatch in vector sizes\n");
	return;
    }

    for (i = 0; i < dest.d; i++) {
	Vn(dest, i) = Vn(src, i);
    }
}

/* Free a vector */
void vec_free(vec_t v) {
    free(v.p);
}
