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

/* vector.h */
/* Routines for dealing with vectors */

#ifndef __vector_h__
#define __vector_h__

#ifdef __cplusplus
extern "C" {
#endif

#include <sys/types.h>

#ifdef WIN32
#include "types.h"
#endif

/* 2D vector of ints */
typedef struct {
    int16_t p[2];
} iv2_t;

/* 3D vector of ints */
typedef struct {
    int16_t p[3];
} iv3_t;

/* 3D vector of ints */
typedef struct {
    int32_t p[3];
} i32_v3_t;

/* 3D vector of ints */
typedef struct {
    int16_t p[5];
} iv5_t;

/* Generalized vector of integers */
typedef struct {
    int d;      /* Dimension */
    int16_t *p; /* Coordinates */
} ivec_t;

/* 2D vector of doubles */
typedef struct {
    double p[2];
} v2_t;

/* 3D vector of doubles */
typedef struct {
    double p[3];
} v3_t;

/* 3D vector of floats */
typedef struct {
    float p[3];
} fv3_t;

/* Generalize vector of doubles */
typedef struct {
    u_int16_t d;      /* Dimension */
    double *p;  /* Coordinates */
} vec_t;

#define Vx(v) (v).p[0]
#define Vy(v) (v).p[1]
#define Vz(v) (v).p[2]

#define Vr(v) (v).p[0]
#define Vg(v) (v).p[1]
#define Vb(v) (v).p[2]

#define Vpx(v) (v)->p[0]
#define Vpy(v) (v)->p[1]
#define Vpz(v) (v)->p[2]

#define Vn(v,n) (v).p[(n)]

/* ******************* 2D vectors ******************* */

/* Constructors for 2D vectors */
iv2_t iv2_new(int16_t x, int16_t y);
v2_t v2_new(double x, double y);

/* Routines for adding 2/3D vectors */
iv2_t iv2_add(iv2_t u, iv2_t v);
v2_t v2_add(v2_t u, v2_t v);
iv2_t iv2_sub(iv2_t u, iv2_t v);
v2_t v2_sub(v2_t u, v2_t v);

/* Scale a vector */
v2_t v2_scale(double s, v2_t v);

/* Compute the norm (length squared) of the given vector */
double v2_norm(v2_t v);

/* Return a unit vector in the same direction as v, v != 0 */
v2_t v2_unit(const v2_t v);

/* Compute the mean of a set of vectors */
v2_t v2_mean(int n, v2_t *v);
void v2_covariance(int n, v2_t *v, v2_t mean, double *cov);
    
/* Dot product */
double v2_dotp(v2_t u, v2_t v);

/* Return the angle between two vectors */
double v2_angle(v2_t u, v2_t v);

/* Compute the pair-wise minimum / maximum of two vectors */
v2_t v2_minimum(v2_t u, v2_t v);
v2_t v2_maximum(v2_t u, v2_t v);

/* Compute the centroid of an array of 2D vectors */
v2_t v2_compute_centroid(v2_t *pts, int num_pts);
iv2_t iv2_compute_centroid(iv2_t *pts, int num_pts);

/* ******************* 3D vectors ******************* */

/* Constructors for 3D vectors */

/* Constructors for 2D vectors */
iv3_t iv3_new(int16_t x, int16_t y, int16_t z);
v3_t v3_new(double x, double y, double z);

/* Add/subtract two 3D vectors */
v3_t v3_add(const v3_t u, const v3_t v);
v3_t v3_sub(const v3_t u, const v3_t v);

/* Scale the given vector by the given amount */
v3_t v3_scale(double s, const v3_t v);

/* Return the (squared/unsquared) magnitude of the given vector */
double v3_magsq(const v3_t v);
double v3_mag(const v3_t v);

/* Compute coordinate-wise min, max of two vectors */
v3_t v3_min(const v3_t u, const v3_t v);
v3_t v3_max(const v3_t u, const v3_t v);

/* Return a unit vector in the same direction as v, v != 0 */
v3_t v3_unit(const v3_t v);

/* Scale the vector so that the 3rd coordinate is 1 */
v3_t v3_homogenize(const v3_t v);

/* Return the dot product of two 3D vectors */
double v3_dotp(const v3_t u, const v3_t v);

/* Return the cross product of two 3D vectors */
v3_t v3_cross(const v3_t u, const v3_t v);

/* Compute the mean of a set of vectors */
v3_t v3_mean(int n, const v3_t *v);
v3_t v3_median(int n, const v3_t *v);
double v3_variance_zm(int n, const v3_t *v);
void v3_covariance_zm(int n, const v3_t *v, double *cov);
void v3_svd(int n, const v3_t *v, double *U, double *S, double *VT);

/* Find the vector in v that is furthest from u */
v3_t v3_extremum(int n, const v3_t *v, const v3_t u);
v3_t v3_extremum2(int n, const v3_t *a, const v3_t u, v3_t v);
int v3_extremum_idx(int n, const v3_t *v, const v3_t u);

/* Assumes unit_normal is normalized */
v3_t v3_project(const v3_t v, const v3_t unit_normal);

/* Print out the vector */
void v3_print(const v3_t v);



/* ******************* ND vectors ******************* */

/* Constructor */
vec_t vec_new(int d);

/* Create the new vector and set all coordinates to the value `v' */
vec_t vec_new_set(int d, double val);

/* Add two vectors */
vec_t vec_add(vec_t u, vec_t v);

/* Subtract two vectors */
vec_t vec_sub(vec_t u, vec_t v);

/* Scale the given vector with the given scalar (changing the
 * given vector) */
void vec_scale_inplace(double c, vec_t v);

/* Compute the norm (length squared) of the given vector */
double vec_norm(vec_t v);

/* Copy one vector into another */
void vec_copy(vec_t dest, vec_t src);

/* Free a vector */
void vec_free(vec_t v);

#ifdef __cplusplus
}
#endif

#endif /* __vector_h__ */
