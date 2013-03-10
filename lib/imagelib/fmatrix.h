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

/* fmatrix.h */
/* Routines for estimating the fundamental matrix of a pair of images */

#ifndef __fmatrix_h__
#define __fmatrix_h__

#ifdef __cplusplus
extern "C" {
#endif

// #include "dmap.h"
#include "image.h"
#include "vector.h"

/* Compute the epipoles of an F-matrix */
void fmatrix_compute_epipoles(double *F, double *e1, double *e2);

/* Compute the distance from l to the epipolar line of r under F */
double fmatrix_compute_residual(double *F, v3_t r, v3_t l);

/* Use RANSAC to estimate an F-matrix */
// void estimate_fmatrix_ransac(img_t *img, img_dmap_t *map, 
//                              int num_trials, double threshold, 
//                              int essential, double *F);

/* Use RANSAC to estimate an F-matrix */
int estimate_fmatrix_ransac_matches(int num_pts, v3_t *a_pts, v3_t *b_pts, 
                                    int num_trials, double threshold, 
                                    double success_ratio, 
                                    int essential, double *F);

/* Use linear least-squares to estimate the fundamantal matrix.  The
 * F-matrix is returned in Fout, and the two epipoles in e1 and e2 */
int estimate_fmatrix_linear(int num_pts, v3_t *r_pts, v3_t *l_pts, 
                            int essential, 
                            double *Fout, double *e1, double *e2);


/* Estimate the essential matrix from an F-matrix, assuming 
 * same focal lengths */
void estimate_essential_same_focal_lengths(double *F, double *alpha, double *E);

/* Estimate the essential matrix from an F-matrix, assuming 
 * different focal lengths */
void estimate_essential_different_focal_lengths(double *F, 
						double *calib1, double *calib2,
						double *E);
    
/* F is the fundamental matrix between image B and image A.  The two
 * homographies that rectify image A and image B are output. */
void fmatrix_rectify_images(img_t *a, img_t *b, double *F, double *H_a, double *H_b);

/* Find the closest rank 2 matrix to the given 3x3 matrix */
int closest_rank2_matrix(double *Fin, double *Fout, double *U, double *VT);

/* Refine an F-matrix estimate using LM */
// void refine_fmatrix_nonlinear(img_dmap_t *amap, double *F0, double *Fout);

/* Refine an F-matrix estimate using LM */
void refine_fmatrix_nonlinear_matches(int num_pts, v3_t *r_pts, v3_t *l_pts, 
				      double *F0, double *Fout);

/* Compute an F-matrix from two sets of camera parameters */
void fmatrix_from_parameters(double *i0, double *R0, double *t0, double *i1, double *R1, double *t1, double *F);

#ifdef __cplusplus
}
#endif

#endif /* __fmatrix_h__ */
