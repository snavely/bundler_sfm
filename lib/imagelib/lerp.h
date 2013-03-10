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

/* lerp.h */
/* Interpolation routines */

#ifndef __lerp_h__
#define __lerp_h__

#define LERP(x0, x1, f0, f1, f2, f3) ((1.0 - (x1)) * ((1.0 - (x0)) * (f0) + (x0) * (f1)) + (x1) * ((1.0 - (x0)) * (f2) + (x0) * (f3)))

/* Interpolate 0 <= x <= 1 between the given function values at points
 * on the [0,1]^n grid.
 *
 * x: an array of length n
 * f: an array of length 2^n
 *
 * To see how this works, see 
 *    http://osl.iu.edu/~tveldhui/papers/MAScThesis/node33.html */
double lerp(int n, double *x, double *f);

/* Special case lerp for 2 dimensions */
double lerp2(double x, double y, int *f);

/* Bilinear interpolation on a 2D function */
double func_lerp(int w, int h, double *f, double x, double y);
double func_lerpf(int w, int h, float *f, double x, double y);

/* Interpolate 0 <= x <= 1 between the given function values at points
 * on the [0,1]^n grid.  This generalized bilinear interpolation */
double nlerp(int n, double *xin, double *f);

double nlerp2(double *xin, double *f);

/* Apply bicubic interpolation to a 4x4 patch around the coordinate
 * (x,y) */
double bicubic_interpolate_2D(double x, double y, double *f);
double bicuberp(double x, double y, double *f);

#endif /* __lerp_h__ */
