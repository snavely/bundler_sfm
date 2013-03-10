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

/* poly3.h */
/* 3-variable polynomial functions */

#ifndef __POLY3_H__
#define __POLY3_H__

typedef enum {
	POLY3_UNIT,
	POLY3_X, POLY3_Y, POLY3_Z, 
	POLY3_XY, POLY3_XZ, POLY3_YZ,
	POLY3_X2, POLY3_Y2, POLY3_Z2,
	POLY3_X2Y, POLY3_X2Z, POLY3_XY2, POLY3_Y2Z, POLY3_XZ2, POLY3_YZ2, POLY3_XYZ,
	POLY3_X3, POLY3_Y3, POLY3_Z3,
	NUM_POLY3_COEFFS
} poly3_coeff_t;

typedef struct {
	double v[NUM_POLY3_COEFFS];
} poly3_t;

poly3_t poly3_new(double x, double y, double z, double unit);
poly3_t poly3_add(poly3_t a, poly3_t b);
poly3_t poly3_add3(poly3_t a, poly3_t b, poly3_t c);
poly3_t poly3_sub(poly3_t a, poly3_t b);
poly3_t poly3_mult(poly3_t a, poly3_t b);
poly3_t poly3_mult11(poly3_t a, poly3_t b);
poly3_t poly3_mult21(poly3_t a, poly3_t b);
poly3_t poly3_scale(poly3_t a, double scale);
double poly3_eval(poly3_t a, double x, double y, double z);
double poly3_get(poly3_t a, poly3_coeff_t idx);

#endif /* __POLY3_H__ */
