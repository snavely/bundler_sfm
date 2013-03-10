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

/* fit.h */
/* Routines for fitting various curves to data */

#ifndef __fit_h__
#define __fit_h__

#ifdef __cplusplus
extern "C" {
#endif

#include "vector.h"

/* Compute the distance between a plane and a point */
double plane_point_distance(double *plane, v3_t pt);
double line_point_distance(double *line, v2_t pt);

/* Project a point onto a plane */
v3_t project_point_onto_plane(double *plane, v3_t pt);
    
/* Fit a line to a set of 2D points */
void fit_2D_line(int num_pts, v2_t *pts, double *params);

/* Fit a line to a set of 2D points */
double fit_2D_line_orthogonal_regression(int num_pts, v2_t *pts, 
					 double *params);

/* Fit a line using orthogonal regression and RANSAC */
double fit_2D_line_ortreg_ransac(int num_pts, v2_t *pts, 
                                 int num_ransac_rounds, 
                                 double ransac_threshold, 
                                 int *num_inliers_out, double *params);

/* Fit a plane to a set of 3D points */
void fit_3D_plane(int num_pts, v3_t *pts, double *params);

/* Fit a plane (which passes through the origin) to a set of 3D points */
void fit_3D_plane_through_origin(int num_pts, v3_t *pts, double *params);
void fit_3D_plane_through_origin2(int num_pts, v3_t *pts, double *params);

/* Fit a plane using orthogonal regression */
double fit_3D_plane_orthogonal_regression(int num_pts, v3_t *pts, 
					  double *params);

/* Fit a plane using orthogonal regression and RANSAC */
double fit_3D_plane_ortreg_ransac(int num_pts, v3_t *pts, 
				  int num_ransac_rounds, 
				  double ransac_threshold, 
				  int *num_inliers_out, double *params);
    
#ifdef __cplusplus
}
#endif

#endif /* __fit_h__ */
